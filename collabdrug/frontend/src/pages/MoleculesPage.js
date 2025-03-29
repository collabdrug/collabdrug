import React from 'react';
import { Container, Typography, Box } from '@mui/material';
import MoleculeUpload from '../components/MoleculeUpload';
import MoleculeList from '../components/MoleculeList';

const MoleculesPage = () => {
  return (
    <Container maxWidth="lg">
      <Box sx={{ py: 4 }}>
        <Typography variant="h4" component="h1" gutterBottom>
          Molecular Data Management
        </Typography>
        
        <MoleculeUpload />
        <MoleculeList />
      </Box>
    </Container>
  );
};

export default MoleculesPage; 