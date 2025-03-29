from pydantic import BaseModel
from datetime import datetime
from typing import List, Optional

class MoleculeBase(BaseModel):
    smiles: str
    name: str
    molecular_weight: float

class MoleculeCreate(MoleculeBase):
    pass

class MoleculePropertyBase(BaseModel):
    property_name: str
    property_value: float

class MoleculePropertyCreate(MoleculePropertyBase):
    pass

class MoleculeProperty(MoleculePropertyBase):
    id: int
    molecule_id: int
    created_at: datetime

    class Config:
        from_attributes = True

class Molecule(MoleculeBase):
    id: int
    created_at: datetime
    updated_at: datetime
    owner_address: Optional[str] = None
    properties: List[MoleculeProperty] = []

    class Config:
        from_attributes = True

class MoleculeResponse(Molecule):
    pass 