import pubchempy as pcp
import time
from pathlib import Path
from urllib.parse import quote

from IPython.display import Markdown, Image
import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import MolsToGridImage


def get_similarity_pubchem(smiles:str, threshold =75, n_records =12):
    """
    Fetch similar molecules from PubChem based on a given SMILES string.
    Args:
        smiles (str): The SMILES representation of the molecule.
        threshold (int, optional): Similarity threshold (0-100). Defaults to 75.
        n_records (int, optional): Number of similar molecules to fetch. Defaults to 12.
    Returns:
    str
    The job key"""
    escaped_smiles = quote(smiles).replace("/", ".")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{escaped_smiles}/JSON?Threshold={threshold}&MaxRecords={n_records}"
    response = requests.get(url)
    response.raise_for_status()
    key = response.json()["Waiting"]["ListKey"]
    return key

def download_pubchem_cids(key:str, attempts = 30):
    """
    Download CIDs from PubChem using a job key.
    Args:
        key (str): The job key obtained from get_similarity_pubchem.
        attempts (int, optional): Number of attempts to check job status. Defaults to 30.
    Returns:
    list
    A list of CIDs."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{key}/cids/JSON"
    print(f"Querying for job {key} at URL {url}...", end="")
    while attempts:
        r = requests.get(url)
        r.raise_for_status()
        response = r.json()
        if "IdentifierList" in response:
            cids = response["IdentifierList"]["CID"]
            break
        attempts -= 1
        print(".", end="")
        time.sleep(10)
    else:
        raise ValueError(f"Could not find matches for job key: {key}")
    return cids

def cid_to_smiles(cids):
    """
    Convert a list of CIDs to their corresponding SMILES strings.
    Args:
        cids (list): List of PubChem CIDs.
    Returns:
    list
    A list of SMILES strings."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(map(str, cids))}/property/CanonicalSMILES,IsomericSMILES/JSON"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    
    if "PropertyTable" not in data or "Properties" not in data["PropertyTable"]:
        raise ValueError("Unexpected response format from PubChem API")
        
    smiles_list = []
    for item in data["PropertyTable"]["Properties"]:
        # Try CanonicalSMILES first, fall back to IsomericSMILES if not available
        if "CanonicalSMILES" in item:
            smiles_list.append(item["CanonicalSMILES"])
        elif "IsomericSMILES" in item:
            smiles_list.append(item["IsomericSMILES"])
        else:
            print(f"Warning: No SMILES found for CID {item.get('CID', 'unknown')}")
    
    return smiles_list
    

if __name__ == "__main__":
    test_smiles = "CCO" # Aspirin
    key = get_similarity_pubchem(test_smiles, threshold=80, n_records=10)
    cids = download_pubchem_cids(key)
    smiles_list = cid_to_smiles(cids)
    print("Similar SMILES from PubChem:")
    for smi in smiles_list:
        print(smi)




