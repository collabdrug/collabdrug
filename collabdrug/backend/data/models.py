from sqlalchemy import Column, Integer, String, Float, DateTime, ForeignKey
from sqlalchemy.orm import relationship
from datetime import datetime
from .database import Base

class Molecule(Base):
    __tablename__ = "molecules"

    id = Column(Integer, primary_key=True, index=True)
    smiles = Column(String, unique=True, index=True)
    name = Column(String)
    molecular_weight = Column(Float)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    owner_address = Column(String)  # Solana wallet address

    properties = relationship("MoleculeProperty", back_populates="molecule")
    predictions = relationship("MoleculePrediction", back_populates="molecule")

class MoleculeProperty(Base):
    __tablename__ = "molecule_properties"

    id = Column(Integer, primary_key=True, index=True)
    molecule_id = Column(Integer, ForeignKey("molecules.id"))
    property_name = Column(String)
    property_value = Column(Float)
    created_at = Column(DateTime, default=datetime.utcnow)

    molecule = relationship("Molecule", back_populates="properties")

class MoleculePrediction(Base):
    __tablename__ = "molecule_predictions"

    id = Column(Integer, primary_key=True, index=True)
    molecule_id = Column(Integer, ForeignKey("molecules.id"))
    prediction_type = Column(String)  # e.g., "drug_target_interaction", "toxicity"
    prediction_value = Column(Float)
    confidence_score = Column(Float)
    created_at = Column(DateTime, default=datetime.utcnow)

    molecule = relationship("Molecule", back_populates="predictions") 