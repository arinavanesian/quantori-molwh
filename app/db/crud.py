from ..models.Users import Users,MoleculeWarehouse
from sqlalchemy.orm import Session
from sqlalchemy import select
from typing import List, Dict, Optional, Tuple, TypeVar

import heapq
import pickle
from ..utility import get_molecules_paginated
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, MACCSkeys

from ..modelbody import (MoleculeRequest, MoleculeResponse,
                         MoleculeListResponse, PaginatedResponse, 
                         SimilarityResponse)

T = TypeVar("T", bound=MoleculeWarehouse)


def get_molecule_by_smiles(
          db_session: Session,
          smiles:str)->Optional[T]:
     return db_session.query(MoleculeWarehouse).filter(MoleculeWarehouse.smiles == smiles).first()

def get_molecule_by_uuid(
          db_session: Session,
          uuid:str)->Optional[T]:
     return db_session.query(MoleculeWarehouse).filter(MoleculeWarehouse.uuid == uuid).first()

def get_molecule_by_iupac(db_session:Session, iupac_name:str)->Optional[T]:
     return db_session.query(MoleculeWarehouse).filter(MoleculeWarehouse.iupac_name==iupac_name).first()

def get_Mols(db_session:Session,skip:int = 0, limit:int=0)->List[Chem.Mol]:
     result = []
     molecules = list_molecules(db_session,skip = skip, limit=limit)

     for mol_db in molecules:
          Mol = pickle.loads(mol_db) if mol_db else None
          result.append((mol_db.uuid,mol_db.smiles, mol_db.iupac_name, Mol))
     return result

def get_Mols_all(db_session:Session):
     return get_Mols(db_session, skip = 0, limit=get_dataCount(db_session))

# TODO: make sure to use as a dependency or a decorator 
def get_cannonical(smiles:str)->str:
     """Returns the canonical SMILES string for a given SMILES input
     making sure that stereoisomers are identifiable"""
     mol = Chem.MolFromSmiles(smiles)
     can_smiles = Chem.MolToSmiles(mol, canonical=True) if mol else None
     return can_smiles

def get_dataCount(db_session:Session)->int:
     count = db_session.query(MoleculeWarehouse).count()
     return count

# TODO: return a tuple of list and it's length 
def list_molecules(
          db_session:Session,
          skip:int = 0,
          limit:int = 50)->List[T]:
     list_molecules = db_session.query(MoleculeWarehouse).offset(skip).limit(limit).all()
     return list_molecules
          
def get_molecules_paginated(db_session: Session, skip: int, limit: int) -> dict:
        """Get paginated molecules with metadata."""
        query = db_session.query(MoleculeWarehouse)
        total_count = query.count()
        molecules = query.offset(skip).limit(limit).all()
        page = (skip // limit) + 1 if limit > 0 else 1
        total_pages = (total_count + limit - 1) // limit if limit > 0 else 1
        return {
            "molecules": molecules,
            "total_count": total_count,
            "page": page,
            "per_page": limit,
            "total_pages": total_pages,
            "has_next": skip + limit < total_count,
            "has_prev": skip > 0}

def get_molecules_paginated(
    db_session: Session, 
    skip: int = 0, 
    limit: int = 100
) -> Dict:
    """
    Get paginated molecules with metadata.
    
    Returns raw database models - this is the CORRECT approach for CRUD layer.
    """
    try:
        query = db_session.query(MoleculeWarehouse)
        total_count = query.count()
        molecules = query.offset(skip).limit(limit).all()
        
        page = (skip // limit) + 1 if limit > 0 else 1
        total_pages = (total_count + limit - 1) // limit if limit > 0 else 1
        
        return {
            "molecules": molecules,
            "total_count": total_count,
            "page": page,
            "per_page": limit,
            "total_pages": total_pages,
            "has_next": skip + limit < total_count,
            "has_prev": skip > 0 
        }
        
    except Exception as e:
        print(f"Error in get_molecules_paginated: {e}")
        raise

#TODO: passing the molreq
def create_molecule(
        db_session: Session,
        uuid: str,
        iupac_name: Optional[str],
        mol: Optional[Chem.Mol],
)->MoleculeWarehouse:
     smiles = Chem.MolToSmiles(mol)
    #  mol_bin = pickle.dumps(mol)
     
     db_molecule = MoleculeWarehouse(
        uuid=uuid,
        smiles=smiles,
        iupac_name=iupac_name,
        # mol_binary=mol_bin
    )
     db_session.add(db_molecule)
     db_session.commit()
     db_session.refresh(db_molecule)
     return db_molecule
    
def create_molecule_req(
        db_session: Session,
        mol_res:MoleculeResponse,
)->MoleculeWarehouse:     
     db_molecule = MoleculeWarehouse(
        uuid=mol_res.uuid,
        smiles=mol_res.smiles,
        iupac_name=mol_res.iupac_name,
        # mol_binary=mol_bin
    )
     db_session.add(db_molecule)
     db_session.commit()
     db_session.refresh(db_molecule)
     return db_molecule

def add_molecule(db_session: Session, molecule: MoleculeWarehouse):
    db_session.add(molecule)
    db_session.commit()
    db_session.refresh(molecule)
    return molecule 

def update_molecule(
    db_session: Session,
    uuid: str,
    smiles: str,
    iupac_name: Optional[str],
    # mol: Chem.Mol
) -> MoleculeWarehouse:
    """Update an existing molecule"""
#     db_molecule = get_molecule_by_uuid(db, uuid)
    db_molecule = get_molecule_by_smiles(db_session=db_session, smiles=smiles)
    if not db_molecule:
        return None
    
    db_molecule.smiles = smiles
    db_molecule.iupac_name = iupac_name
    # db_molecule.mol_binary = pickle.dumps(mol) if mol else None
    
    db_session.commit()
    db_session.refresh(db_molecule)
    return db_molecule

def delete_molecule(db_session: Session, uuid:str)->bool:
     db_molecule = get_molecule_by_uuid(db_session, uuid)
     if db_molecule:
         db_session.delete(db_molecule)
         db_session.commit()
         return True
     return False

def get_all(db_session: Session)->List[T]:
     count = get_dataCount(db_session)
     return list_molecules(db_session, skip=0, limit=count)

def substructure_search(db_session: Session, query_smiles:str)->MoleculeListResponse:
     query_mol = Chem.MolFromSmiles(query_smiles)
     molecules = get_Mols_all(db_session=db_session)
     matches =[]
     for uuid, smiles, iupac_name, mol in molecules:
          if mol and mol.HasSubstructMatch(query_mol):
               matches.append(MoleculeResponse(
                    uuid=uuid,
                    smiles=smiles,
                    iupac_name=iupac_name,
                    # mol_bin = pickle.dumps(mol)
               ))
     return MoleculeListResponse(molecules=matches, count=len(matches))

# TODO: return also similarity scores
def similarity_search(db_session: Session, 
                      query_smiles:str,
                        radius:int =  2, 
                      n_bits:int = 2048, 
                      threshold:float=0.7,
                      top_k:int =50)-> List[SimilarityResponse]:
     
     """
     Perform a similarity search using Morgan fingerprints as
     bit vectors and Tanimoto similarity."""
     fpg = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
     query_mol = Chem.MolFromSmiles(query_smiles)
     if not query_mol:
          return MoleculeListResponse(molecules=[], count=0)
     
     query_circular_bit_fp = fpg.GetFingerprint(query_mol)

     molecules = get_Mols_all(db_session=db_session)
     heap: list[tuple[float, MoleculeResponse]] = []
     for uuid, smiles, iupac_name, mol in molecules:
          if not mol:
               continue
          else:
               target_circular_bit_fp = fpg.GetFingerprint(mol)
               similarity = DataStructs.TanimotoSimilarity(query_circular_bit_fp,
                                                     target_circular_bit_fp)
               if similarity >= threshold:
                    mol_resp =MoleculeResponse(
                         uuid=uuid,
                         smiles=smiles,
                         iupac_name=iupac_name,
                         # mol_bin = pickle.dumps(mol)
                         )
                    if len(heap) < top_k:
                         heapq.heappush(heap, (similarity, mol_resp))
                    else:
                         heapq.heappushpop(heap, (similarity, mol_resp))

     heap.sort(key = lambda x: x[0], reverse=True)
     results = [
          SimilarityResponse(
               molecule=mol_resp, similarity_score=sim)
               for sim, mol_resp in heap
     ]     
     return results