from rdkit import Chem
from rdkit.Chem import MACCSkeys

# === Utility ===

def generate_MACCSkeys(query_smiles: str):
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise ValueError(f"Invalid SMILES: {query_smiles!r}")
    maccs_fp = MACCSkeys.GenMACCSKeys(query_mol)
    return maccs_fp

# def generate_MACCSKeys_list(query_list:List[str]):