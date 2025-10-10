# Whole-slide IMC Workflows

This repository provides a lightweight workflow for generating tissue masks from imaging mass cytometry (IMC) data preprocessing and analysis, including:

**1. Tissue Mask Generation:** creates binary masks using dynamic thresholding methods (Otsu or 2-component GuassianMixture Models) with optional small obkect removal and hole filling

**2. Patch Extraction:** generates image patches for downstream analysis, with flexible marker selection, size, and stride

**3. Preprocessing:** includes customizable standardized preprocessing (normalization, scaling, channel selection) automatically applied if not already performed

## Table of Contents
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Arguments](#arguments)
- [Outputs](#outputs)

## Installation

### Using Conda
'''
git clone https://github.com/sophiamjiali/ws_imc_workflows.git
cd ws_imc_workflows
conda env create -f environment.yml
conda activate ws_imc_workflows_env
'''

### Using pip
'''
python -m pip install -r requirements.txt
'''

## Configuration

## Usage

## Arguments

## Outputs

