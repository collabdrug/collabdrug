from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np

def calculate_molecular_weight(mol: Chem.Mol) -> float:
    """Calculate molecular weight of a molecule."""
    return Descriptors.ExactMolWt(mol)

def calculate_logp(mol: Chem.Mol) -> float:
    """Calculate logP of a molecule."""
    return Descriptors.MolLogP(mol)

def calculate_hbd(mol: Chem.Mol) -> int:
    """Calculate number of hydrogen bond donors."""
    return Descriptors.NumHDonors(mol)

def calculate_hba(mol: Chem.Mol) -> int:
    """Calculate number of hydrogen bond acceptors."""
    return Descriptors.NumHAcceptors(mol)

def calculate_rotatable_bonds(mol: Chem.Mol) -> int:
    """Calculate number of rotatable bonds."""
    return Descriptors.NumRotatableBonds(mol)

def calculate_psa(mol: Chem.Mol) -> float:
    """Calculate polar surface area."""
    return Descriptors.TPSA(mol)

def calculate_all_properties(mol: Chem.Mol) -> dict:
    """Calculate all molecular properties."""
    return {
        "molecular_weight": calculate_molecular_weight(mol),
        "logp": calculate_logp(mol),
        "hbd": calculate_hbd(mol),
        "hba": calculate_hba(mol),
        "rotatable_bonds": calculate_rotatable_bonds(mol),
        "psa": calculate_psa(mol)
    }

def validate_smiles(smiles: str) -> bool:
    """Validate SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def canonicalize_smiles(smiles: str) -> str:
    """Canonicalize SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    return Chem.MolToSmiles(mol, canonical=True) 