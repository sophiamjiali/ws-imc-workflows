# Whole-slide Imaging Mass Cytometry (WS-IMC) Workflows

This repository provides reproducible workflows for imaging mass cytometry (WS-IMC) data preprocessing and analysis, including:

**1. Tissue Mask Generation:** creates binary masks using dynamic thresholding methods (Otsu or GuassianMixture Models) with optional small obkect removal and hole filling

**2. Patch Extraction:** generates image patches for downstream analysis, with flexible marker selection, size, and stride

Features include:
- YAML-based configuration for all input/output paths and parameters
- Shared utilities for loading images, panel files, and metadata
- Quality control visualization of masks overlaid on original images
- Metadata recording per image for reproducibility and downstream analysis
