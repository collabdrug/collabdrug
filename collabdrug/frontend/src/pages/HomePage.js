import React from 'react';
import { Container, Typography, Box, Button, Grid, Paper } from '@mui/material';
import { Link as RouterLink } from 'react-router-dom';
import { Science as ScienceIcon, Storage as StorageIcon, Timeline as TimelineIcon } from '@mui/icons-material';

const HomePage = () => {
  return (
    <Container maxWidth="lg">
      <Box sx={{ py: 8 }}>
        {/* Hero Section */}
        <Box sx={{ textAlign: 'center', mb: 8 }}>
          <Typography variant="h2" component="h1" gutterBottom>
            Welcome to CollabDrug
          </Typography>
          <Typography variant="h5" color="text.secondary" paragraph>
            A decentralized platform for collaborative drug discovery
          </Typography>
          <Button
            variant="contained"
            size="large"
            component={RouterLink}
            to="/molecules"
            sx={{ mt: 2 }}
          >
            Get Started
          </Button>
        </Box>

        {/* Features Section */}
        <Grid container spacing={4}>
          <Grid item xs={12} md={4}>
            <Paper sx={{ p: 3, height: '100%', textAlign: 'center' }}>
              <ScienceIcon sx={{ fontSize: 48, color: 'primary.main', mb: 2 }} />
              <Typography variant="h6" gutterBottom>
                Molecular Data Management
              </Typography>
              <Typography color="text.secondary">
                Upload, manage, and analyze molecular data with ease
              </Typography>
            </Paper>
          </Grid>
          <Grid item xs={12} md={4}>
            <Paper sx={{ p: 3, height: '100%', textAlign: 'center' }}>
              <StorageIcon sx={{ fontSize: 48, color: 'primary.main', mb: 2 }} />
              <Typography variant="h6" gutterBottom>
                Blockchain Integration
              </Typography>
              <Typography color="text.secondary">
                Secure data ownership and transparent contribution tracking
              </Typography>
            </Paper>
          </Grid>
          <Grid item xs={12} md={4}>
            <Paper sx={{ p: 3, height: '100%', textAlign: 'center' }}>
              <TimelineIcon sx={{ fontSize: 48, color: 'primary.main', mb: 2 }} />
              <Typography variant="h6" gutterBottom>
                AI-Powered Predictions
              </Typography>
              <Typography color="text.secondary">
                Leverage advanced AI models for drug discovery
              </Typography>
            </Paper>
          </Grid>
        </Grid>
      </Box>
    </Container>
  );
};

export default HomePage; 