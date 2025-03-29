"""initial

Revision ID: 001
Revises: 
Create Date: 2024-01-01 00:00:00.000000

"""
from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = '001'
down_revision = None
branch_labels = None
depends_on = None

def upgrade() -> None:
    # Create molecules table
    op.create_table(
        'molecules',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('smiles', sa.String(), nullable=False),
        sa.Column('name', sa.String(), nullable=False),
        sa.Column('molecular_weight', sa.Float(), nullable=False),
        sa.Column('created_at', sa.DateTime(), nullable=False),
        sa.Column('updated_at', sa.DateTime(), nullable=False),
        sa.Column('owner_address', sa.String(), nullable=True),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_molecules_id'), 'molecules', ['id'], unique=False)
    op.create_index(op.f('ix_molecules_smiles'), 'molecules', ['smiles'], unique=True)

    # Create molecule_properties table
    op.create_table(
        'molecule_properties',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('molecule_id', sa.Integer(), nullable=False),
        sa.Column('property_name', sa.String(), nullable=False),
        sa.Column('property_value', sa.Float(), nullable=False),
        sa.Column('created_at', sa.DateTime(), nullable=False),
        sa.ForeignKeyConstraint(['molecule_id'], ['molecules.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_molecule_properties_id'), 'molecule_properties', ['id'], unique=False)

    # Create molecule_predictions table
    op.create_table(
        'molecule_predictions',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('molecule_id', sa.Integer(), nullable=False),
        sa.Column('prediction_type', sa.String(), nullable=False),
        sa.Column('prediction_value', sa.Float(), nullable=False),
        sa.Column('confidence_score', sa.Float(), nullable=False),
        sa.Column('created_at', sa.DateTime(), nullable=False),
        sa.ForeignKeyConstraint(['molecule_id'], ['molecules.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_molecule_predictions_id'), 'molecule_predictions', ['id'], unique=False)

def downgrade() -> None:
    op.drop_index(op.f('ix_molecule_predictions_id'), table_name='molecule_predictions')
    op.drop_table('molecule_predictions')
    
    op.drop_index(op.f('ix_molecule_properties_id'), table_name='molecule_properties')
    op.drop_table('molecule_properties')
    
    op.drop_index(op.f('ix_molecules_smiles'), table_name='molecules')
    op.drop_index(op.f('ix_molecules_id'), table_name='molecules')
    op.drop_table('molecules') 