from rdkit import Chem
from rdkit.Chem import MACCSkeys
from sqlalchemy.orm import Session


def generate_MACCSkeys(query_smiles: str):
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise ValueError(f"Invalid SMILES: {query_smiles!r}")
    maccs_fp = MACCSkeys.GenMACCSKeys(query_mol)
    return maccs_fp