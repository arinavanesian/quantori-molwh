import pytest
from fastapi.testclient import TestClient
from main import app
from fastapi import FastAPI

client = TestClient(app)

@pytest.fixture(autouse=True)
def setup():
    """Create fresh store before each test"""
    client.post("/store/create")

def test_health_check():
    """Test health endpoint works"""
    response = client.get("/health")
    assert response.status_code == 200

def test_add_molecule():
    """Test adding a molecule"""
    response = client.post("/add", json={"smiles": "CCO"})
    assert response.status_code == 201
    data = response.json()
    assert data["smiles"] == "CCO"
    assert "uuid" in data

def test_list_empty():
    """Test listing when store is empty"""
    response = client.get("/list")
    assert response.status_code == 200
    data = response.json()
    assert data["count"] == 0

def test_list_with_molecules():
    """Test listing after adding molecules"""
    client.post("/add", json={"smiles": "CCO"})
    response = client.get("/list")
    assert response.status_code == 200
    data = response.json()
    assert data["count"] == 1

def test_search():
    """Test substructure search"""
    client.post("/add", json={"smiles": "CCO"})
    response = client.post("/search", params={"smiles": "CC"})
    assert response.status_code == 200
    matches = response.json()
    assert isinstance(matches, list)


def test_invalid_smiles():
    """Test adding invalid SMILES fails"""
    response = client.post("/add", json={"smiles": "invalid"})
    assert response.status_code == 422

def test_molecule_not_found():
    """Test getting non-existent molecule"""
    response = client.get("/get/fake_name")
    assert response.status_code == 404