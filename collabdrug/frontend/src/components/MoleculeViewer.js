import React, { useEffect, useRef } from 'react';
import { Box, Paper, Typography } from '@mui/material';

const MoleculeViewer = ({ smiles, width = 300, height = 300 }) => {
  const canvasRef = useRef(null);

  useEffect(() => {
    const drawMolecule = async () => {
      if (!smiles || !canvasRef.current) return;

      try {
        // Initialize RDKit
        const rdkit = await window.rdkitModule();
        
        // Create molecule from SMILES
        const mol = rdkit.get_mol(smiles);
        
        // Get canvas context
        const canvas = canvasRef.current;
        const ctx = canvas.getContext('2d');
        
        // Clear canvas
        ctx.clearRect(0, 0, width, height);
        
        // Draw molecule
        const svg = mol.get_svg_with_highlights(JSON.stringify({}));
        
        // Create SVG element
        const svgElement = new DOMParser().parseFromString(svg, 'image/svg+xml').documentElement;
        
        // Convert SVG to canvas
        const img = new Image();
        const svgBlob = new Blob([svg], { type: 'image/svg+xml;charset=utf-8' });
        const url = URL.createObjectURL(svgBlob);
        
        img.onload = () => {
          ctx.drawImage(img, 0, 0, width, height);
          URL.revokeObjectURL(url);
        };
        
        img.src = url;
        
        // Cleanup
        mol.delete();
      } catch (error) {
        console.error('Error drawing molecule:', error);
      }
    };

    drawMolecule();
  }, [smiles, width, height]);

  return (
    <Paper sx={{ p: 2, textAlign: 'center' }}>
      <canvas
        ref={canvasRef}
        width={width}
        height={height}
        style={{ border: '1px solid #ddd', borderRadius: '4px' }}
      />
      <Typography variant="body2" sx={{ mt: 1, wordBreak: 'break-all' }}>
        {smiles}
      </Typography>
    </Paper>
  );
};

export default MoleculeViewer; 