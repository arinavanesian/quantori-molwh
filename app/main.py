from fastapi import FastAPI, HTTPException, Depends, Query, status
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field, field_validator
from typing import Dict, List, Optional, Literal
import bisect
import logging
from contextlib import asynccontextmanager
import tempfile
import os

from sqlalchemy.orm import Session
from sqlalchemy import text 

from .db.engine import engine
from .db.Base import Base
from .db.session import session_scope, get_db
from .db import crud

from .modelbody import (
    MoleculeResponse, 
    MoleculeListResponse, 
    MoleculeRequest,
    PaginatedResponse, 
    SimilarityResponse)

from .models.Users import MoleculeWarehouse
from .auth import routes as auth_routes

from rdkit import Chem
from rdkit.Chem import MACCSkeys
# from rdkit.Chem import Draw
from rdkit import DataStructs
import pubchempy as pcp

import pandas as pd
import uuid

import uvicorn


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize the application state and create the molecule store"""
    logger.info("Creating the MoleculeWH-DB tables ...    ")
    Base.metadata.create_all(bind=engine)
    logger.info("Molecule DB created successfully.")
    yield
    logger.info("Application shutting down...")

__version__ = "1.0.0"
app = FastAPI(
    version=__version__,
    title=f"Molecule Substructure Search Warehouse v{__version__}",
    description="API for storing, searching, and analyzing molecular structures",
    lifespan=lifespan
)

app.include_router(auth_routes.router)


# ===== Dependencies =====

# TODO: call this in the function
def validate_smiles(smiles:str) -> Chem.Mol:
    """Validate a SMILES string and return an RDKit Mol object"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
            raise HTTPException(
                status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
                detail="Invalid SMILES string"
            )
    return mol

def get_MolFromSmiles(smiles:str)->Chem.Mol:
    return validate_smiles(smiles)

def get_iupac_name(smiles_string: str) -> Optional[str]:
    
    """
    Retrieves the IUPAC name for a given SMILES string using PubChem.
    Args:
        smiles_string (str): The SMILES string of the molecule.
    Returns:
        Optional[str]: The IUPAC name of the compound, or None if not found.
    """
    try:
        mol =validate_smiles(smiles_string)     
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        
        compounds = pcp.get_compounds(canonical_smiles, 'smiles')
        if compounds and compounds[0].iupac_name:
            logger.info(f"Found IUPAC name '{compounds[0].iupac_name}' for SMILES {smiles_string}")
            return compounds[0].iupac_name
        return None
        
    except pcp.PubChemPyError as e:
        logger.warning(f"PubChem lookup failed for SMILES {smiles_string}: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error in IUPAC lookup: {e}")
        return None

# TODO: Check
def validate_existing_molecule(
        smiles: str,
        db_session: Session = Depends(get_db))->str:
    """Validate that a molecule exists in the store
    using either its IUPAC name or UUID as the key."""
    checking_molecule = crud.get_molecule_by_smiles(db_session, smiles)
    if not checking_molecule:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Molecule with '{smiles}' not found"
        )
    return checking_molecule



# @app.post("/store/create", status_code=status.HTTP_201_CREATED)
# def create_store():
#     """Create a new empty molecule store"""
#     global molecule_store
#     molecule_store = pd.DataFrame(columns=["uuid", "iupac_name", "smiles", "Mol"]).set_index("uuid")
#     logger.info("New molecule store created")
#     return {"detail": "Molecule store created successfully!"}

# TODO: add dependency validate_existing
@app.post("/add", response_model=MoleculeResponse, status_code=status.HTTP_201_CREATED)
def add_molecule(
    mol_req: MoleculeRequest,
    db_session:Session = Depends(get_db)
):
    """Adds a new molecule to the store """
    try:
        existing = crud.get_molecule_by_smiles(db_session, mol_req.smiles)
        if existing:
            raise HTTPException(
                status_code=status.HTTP_409_CONFLICT,
                detail="Molecule with the same SMILES already exists"
            )
        validated_mol = validate_smiles(mol_req.smiles)
    
        molecule_uuid = str(uuid.uuid4())
        iupac_name = get_iupac_name(mol_req.smiles) 
        added_molecule = crud.create_molecule(
            db_session=db_session,
            uuid=molecule_uuid,
            iupac_name=iupac_name,
            mol=validated_mol
        )
        mol_resp = MoleculeResponse(
            uuid=molecule_uuid,
            smiles=mol_req.smiles,
            iupac_name=iupac_name
            # mol_bin = pickle.dumps(Mol)
        )
        logger.info(f"Added molecule {molecule_uuid} to store")    
        return mol_resp
    except HTTPException:
        raise       
    except Exception as e:
        logger.error(f"Error adding molecule: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to add molecule to store"
        )

@app.get("/get/{smiles}", response_model=MoleculeResponse)
def get_molecule(
    smiles: str,
    db_session: Session = Depends(get_db)):
    """Get a molecule by its key (UUID or IUPAC name)"""
    db_molecule = crud.get_molecule_by_smiles(db_session=db_session, smiles=smiles)

    return MoleculeResponse(
        uuid=db_molecule.uuid,
        smiles=db_molecule.smiles,
        iupac_name=db_molecule.iupac_name
    )


@app.get("/iupac/{smiles}")
def fetch_iupac(smiles: str):
    iupac_name = get_iupac_name(smiles)
    return {"iupac_name": iupac_name, "smiles": smiles}

# TODO:
@app.put("/upload/{iupac_name}", response_model=MoleculeResponse)
def update_molecule(
    smiles: str,
    mol: MoleculeRequest = ...,
    db_session:Session = Depends(get_db)
):
    """Update an existing molecule by  key SMILES"""
    try:
        existing = crud.get_molecule_by_smiles(db_session, smiles)
        if existing:
            if existing.smiles==mol.smiles and existing.iupac_name==mol.iupac_name:
                mol_resp = MoleculeResponse(
                uuid=uuid,
                smiles=smiles,
                iupac_name=mol.iupac_name
                # mol_bin = pickle.dumps(Chem.MolFromSmiles(existing.smiles))
            )
            
        logger.info(f"Updated molecule {uuid}")
        return mol_resp
        
    except Exception as e:
        logger.error(f"Error updating molecule {uuid}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to update molecule"
        )
    
# TODO:can change to a smiles or iupac_name
# TODO: admin only
@app.delete("/delete/{mol_uuid}", status_code=status.HTTP_204_NO_CONTENT)
def delete_molecule(
    mol_uuid: str = Depends(validate_existing_molecule),
    db_session: Session = Depends(get_db)
):
    """Delete a molecule by psasing a UUID"""
    deleted = crud.delete_molecule(db_session, mol_uuid)
    if not deleted:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to delete molecule with {mol_uuid}"
        )
    logger.info(f"Deleted molecule {mol_uuid}")

@app.get("/list_paged", response_model=PaginatedResponse)
def list_molecules_paginated(
          db_session:Session = Depends(get_db),
          skip:int = 0,
          limit:int = 50):
        try:
            result = crud.get_molecules_paginated(db_session, skip, limit)
            molecules_list = [
                MoleculeResponse(
                    uuid=mol.uuid,
                    smiles=mol.smiles,
                    iupac_name=mol.iupac_name,
                    # mol_bin = pickle.dumps(mol)
                ) 
                for mol in result["molecules"]
            ]
            
            return PaginatedResponse(
                total_count=result["total_count"],
                molecules=molecules_list,
                page=result["page"],
                per_page=result["per_page"],
                total_pages=result["total_pages"],
                has_next=result["has_next"],
                has_prev=result["has_prev"])
        except Exception as e:
            logger.error(f"Error in get_molecules_paginated: {e}")
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Failed to retrieve paginated list of molecules"
            )
      
     
@app.get("/molecules/", response_model=MoleculeListResponse)
def list_molecules(
    skip: int = Query(0, ge=0, description="Number of records to skip"),
    limit: int = Query(100, ge=1, le=1000, description="Number of records to return"),
    db_session: Session = Depends(get_db)
)->MoleculeListResponse:
    """List molecules with pagination and count"""
    try:
        result = crud.list_molecules(db_session, skip, limit)
        molecules = [
            MoleculeResponse(
                uuid=mol.uuid,
                smiles=mol.smiles,
                iupac_name=mol.iupac_name,
                # mol_bin = pickle.dumps(mol)
            ) 
            for mol in result
        ]
        total_count = crud.get_dataCount(db_session)
        return MoleculeListResponse(molecules=molecules, count=total_count)
    except Exception as e:
        logger.error(f"Error in list_molecules: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to retrieve list of molecules"
        )
        

@app.post("/search/{query_smiles}", response_model=MoleculeListResponse)
def substructure_search(
    query_smiles: str,
    db_session: Session = Depends(get_db),
):
    """Search for molecules containing the query substructure"""
    query_mol = validate_smiles(query_smiles)
    mol_list_resp = crud.substructure_search(db_session, query_smiles)
   
    logger.info(f"Substructure search found {mol_list_resp.count} matches for query {query_smiles}")
    return mol_list_resp


#TODO: add fingerprint choice (morgan, rdkit, maccs)

@app.post("/search/similarity", response_model=SimilarityResponse)
def similarity_search(
    query_smiles: str,
    db_session: Session = Depends(get_db),
    threshold: float = 0.0,
    limit:int = 10,
    top_k:int =50
)-> List[SimilarityResponse]:
    """Find the most similar molecules based on Tanimoto similarity"""
   
    
    query_mol = validate_smiles(query_smiles)
    # TODO: remove
    matches = crud.similarity_search(
        db_session=db_session,
        query_smiles=query_smiles,
        threshold=threshold,

    )
    if not matches:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"No similar molecules found with similarity >={threshold}"
        )
    
    return matches




@app.get("/health")
def health_check(db_session: Session = Depends(get_db)):
    """Health check endpoint"""
    try:
        db_session.execute(text("SELECT 1"))  
        total_count = db_session.query(MoleculeWarehouse).count()
        return {
            "status": "healthy", 
            "molecule_count": total_count,
            "database": "connected"
        }
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail=f"Database connection failed: {str(e)}"
        )

    


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
