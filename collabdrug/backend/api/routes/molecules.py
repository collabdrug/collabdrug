from fastapi import APIRouter, Depends, HTTPException, UploadFile, File
from sqlalchemy.orm import Session
from typing import List
import rdkit
from rdkit import Chem
from ..schemas.molecule import MoleculeCreate, MoleculeResponse
from ...data.database import get_db
from ...data.models import Molecule, MoleculeProperty
from ...services.molecule_service import calculate_molecular_weight

router = APIRouter()

@router.post("/upload", response_model=MoleculeResponse)
async def upload_molecule(
    file: UploadFile = File(...),
    db: Session = Depends(get_db)
):
    # Read SMILES from file
    content = await file.read()
    smiles = content.decode().strip()
    
    # Validate SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    # Calculate molecular weight
    mw = calculate_molecular_weight(mol)
    
    # Create molecule in database
    db_molecule = Molecule(
        smiles=smiles,
        name=file.filename,
        molecular_weight=mw
    )
    
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)
    
    return db_molecule

@router.get("/{molecule_id}", response_model=MoleculeResponse)
def get_molecule(molecule_id: int, db: Session = Depends(get_db)):
    molecule = db.query(Molecule).filter(Molecule.id == molecule_id).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule

@router.get("/", response_model=List[MoleculeResponse])
def list_molecules(
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db)
):
    molecules = db.query(Molecule).offset(skip).limit(limit).all()
    return molecules

@router.delete("/{molecule_id}")
def delete_molecule(molecule_id: int, db: Session = Depends(get_db)):
    molecule = db.query(Molecule).filter(Molecule.id == molecule_id).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    
    db.delete(molecule)
    db.commit()
    return {"message": "Molecule deleted successfully"} 