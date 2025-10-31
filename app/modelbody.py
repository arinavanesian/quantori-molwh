from pydantic import BaseModel, Field, field_validator
from typing import Dict, List, Optional
from rdkit import Chem

class MoleculeResponse(BaseModel):
    uuid: str
    smiles: str
    iupac_name: Optional[str] = None
    # mol_binary: Optional[bytes] = Field(None, aliads = "molbin", description="Binary representation of the molecule")
    molecular_weight: Optional[float] = Field(None,alias = "molwt", description="Calculated Molecular Weight")
    log_p: Optional[float] = Field(None, alias ="logp",  description="Calculated Octanol-Water Partition Coefficient")
    tpsa: Optional[float] = Field(None, alias = "tpsa", description="Calculated Topological Polar Surface Area")
    
class MoleculeListResponse(BaseModel):
    molecules: List[MoleculeResponse]
    count: int

class MoleculeRequest(BaseModel):
    smiles: str = Field(..., description="SMILES representation of the molecule", min_length=1)
    iupac_name: Optional[str] = Field(None, description="IUPAC name of the molecule")

    @field_validator('smiles')
    def validate_smiles(cls, v):
        if Chem.MolFromSmiles(v) is None:
            raise ValueError('Invalid SMILES string')
        return v

class PaginatedResponse(BaseModel):
    total_count: int
    molecules: List[MoleculeResponse]
    page:int
    per_page:int
    total_pages:int
    has_next: bool
    has_prev: bool

    

class SimilarityResponse(BaseModel):
    molecule: MoleculeResponse
    similarity_score: float