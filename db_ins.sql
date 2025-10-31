-- Insert users only if they don't exist
INSERT INTO users (email, role, is_active)
SELECT 'admin@molwarehouse.com', 'admin', 1
WHERE NOT EXISTS (SELECT 1 FROM users WHERE email = 'admin@molwarehouse.com');

INSERT INTO users (email, role, is_active)
SELECT 'user1@example.com', 'user', 1
WHERE NOT EXISTS (SELECT 1 FROM users WHERE email = 'user1@example.com');

INSERT INTO users (email, role, is_active)
SELECT 'user2@example.com', 'user', 1
WHERE NOT EXISTS (SELECT 1 FROM users WHERE email = 'user2@example.com');

-- Insert molecules only if they don't exist
INSERT INTO molecule_store (uuid, smiles, iupac_name)
SELECT 'a1b2c3d4-e5f6-7890-abcd-ef1234567890', 'CCO', 'ethanol'
WHERE NOT EXISTS (SELECT 1 FROM molecule_store WHERE smiles = 'CCO');

INSERT INTO molecule_store (uuid, smiles, iupac_name)
SELECT 'b2c3d4e5-f6g7-8901-bcde-f23456789012', 'CC(=O)O', 'acetic acid'
WHERE NOT EXISTS (SELECT 1 FROM molecule_store WHERE smiles = 'CC(=O)O');

INSERT INTO molecule_store (uuid, smiles, iupac_name)
SELECT 'c3d4e5f6-g7h8-9012-cdef-345678901234', 'C1CCCCC1', 'cyclohexane'
WHERE NOT EXISTS (SELECT 1 FROM molecule_store WHERE smiles = 'C1CCCCC1');

INSERT INTO molecule_store (uuid, smiles, iupac_name)
SELECT 'd4e5f6g7-h8i9-0123-def0-456789012345', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'caffeine'
WHERE NOT EXISTS (SELECT 1 FROM molecule_store WHERE smiles = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C');