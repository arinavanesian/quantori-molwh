-- Create tables
-- CREATE TABLE users (
--     id INTEGER PRIMARY KEY AUTOINCREMENT,
--     email VARCHAR NOT NULL UNIQUE,
--     role VARCHAR CHECK (role IN ('admin', 'user')) DEFAULT 'user' NOT NULL,
--     is_active BOOLEAN DEFAULT TRUE NOT NULL
-- );

-- CREATE TABLE molecule_store (
--     uuid VARCHAR PRIMARY KEY,
--     smiles VARCHAR NOT NULL UNIQUE,
--     iupac_name VARCHAR
-- );

-- CREATE TABLE query (
--     query_id INTEGER PRIMARY KEY AUTOINCREMENT,
--     user_id INTEGER NOT NULL,
--     molecule_smiles VARCHAR NOT NULL,
--     created_at DATETIME DEFAULT CURRENT_TIMESTAMP NOT NULL,
--     FOREIGN KEY (user_id) REFERENCES users (id),
--     FOREIGN KEY (molecule_smiles) REFERENCES molecule_store (smiles)
-- );

-- CREATE TABLE user_db (
--     user_id INTEGER PRIMARY KEY,
--     hashed_password VARCHAR NOT NULL,
--     FOREIGN KEY (user_id) REFERENCES users (id)
-- );

-- Create indexes
-- CREATE INDEX ix_users_id ON users (id);
-- CREATE INDEX ix_users_email ON users (email);
-- CREATE INDEX ix_molecule_store_uuid ON molecule_store (uuid);
-- CREATE INDEX ix_molecule_store_smiles ON molecule_store (smiles);
-- CREATE INDEX ix_query_query_id ON query (query_id);
-- CREATE INDEX ix_query_user_id ON query (user_id);
-- CREATE INDEX ix_query_molecule_smiles ON query (molecule_smiles);
-- CREATE INDEX ix_user_db_user_id ON user_db (user_id);

-- Insert sample data
-- Users
-- INSERT INTO users (email, role, is_active) VALUES 
-- ('admin@molwarehouse.com', 'admin', TRUE),
-- ('user1@example.com', 'user', TRUE),
-- ('user2@example.com', 'user', TRUE);

-- -- Molecules
-- INSERT INTO molecule_store (uuid, smiles, iupac_name) VALUES 
-- ('CCO', 'ethanol'),
-- ('CC(=O)O', 'acetic acid'),
-- ('C1CCCCC1', 'cyclohexane'),
-- ('CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'caffeine');

-- Users (with INSERT OR IGNORE to handle duplicates)
INSERT OR IGNORE INTO users (email, role, is_active) VALUES 
('admin@molwarehouse.com', 'admin', 1),
('user1@example.com', 'user', 1),
('user2@example.com', 'user', 1);

-- Molecules (let SQLite auto-generate UUIDs, only provide smiles and iupac_name)
INSERT OR IGNORE INTO molecule_store (smiles, iupac_name) VALUES 
('CCO', 'ethanol'),
('CC(=O)O', 'acetic acid'),
('C1CCCCC1', 'cyclohexane'),
('CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'caffeine');

-- UserDB (passwords are hashed 'password123')
-- INSERT INTO user_db (user_id, hashed_password) VALUES 
-- (1, '$2b$12$EixZaYVK1fsbw1ZfbX3OXePaWxn96p36WQoeG6Lruj3vjPGga31lW'),
-- (2, '$2b$12$EixZaYVK1fsbw1ZfbX3OXePaWxn96p36WQoeG6Lruj3vjPGga31lW'),
-- (3, '$2b$12$EixZaYVK1fsbw1ZfbX3OXePaWxn96p36WQoeG6Lruj3vjPGga31lW');

-- -- Query history
-- INSERT INTO query (user_id, molecule_smiles, created_at) VALUES 
-- (1, 'CCO', '2024-01-15 10:30:00'),
-- (2, 'CC(=O)O', '2024-01-15 11:15:00'),
-- (1, 'C1CCCCC1', '2024-01-15 14:20:00'),
-- (3, 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', '2024-01-16 09:45:00');
