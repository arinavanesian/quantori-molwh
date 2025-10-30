import pandas as pd
import pubchempy as pcp

# Get 1000 molecules from PubChem
import pubchempy as pcp

compounds = pcp.get_compounds('aspirin', 'name', listkey_count=1000)
df = pd.DataFrame([{
    'smiles': c.canonical_smiles,
    'iupac_name': c.iupac_name
} for c in compounds])
df.to_csv('molecules.csv', index=False)