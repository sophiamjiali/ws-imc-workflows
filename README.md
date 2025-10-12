# Whole-slide IMC Workflows

This repository provides a lightweight workflow for generating tissue masks from imaging mass cytometry (IMC) data preprocessing and analysis, including:

**1. Tissue Mask Generation:** creates binary masks using dynamic thresholding methods (Otsu or 2-component GuassianMixture Models) with optional preprocessing, mask cleaning (small object removal and/or hole filling), mask generation metadata, and quality control plots.

**2. Patch Extraction:** extracts patches from images as a `.zarr` with optional marker selection and custom patch size, stride, and minimum tissue coverage per patch.

**3. Preprocessing:** includes customizable standardized preprocessing (normalization, scaling, channel selection) automatically applied if not already performed.

## Table of Contents
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)

## Repository Structure
```
ws_imc_workflows/
├── config/                        # YAML configuration files defining workflow parameters
│   ├── template_config.yaml       # - blank template for configurations 
│   └── default.yaml               # - example configurations with default values
├── docs/                          # Documentation regarding the workflow
│   └── METHODS.md                 # - methods for the workflows
├── src/                           # Executable Python scripts for each workflow:
│   ├── extract_patches.ipynb      # - exporatory patch extraction (Jupyter Notebook)
│   ├── extract_patches.py         # - patch extraction workflow
│   ├── generate_masks.ipynb       # - exploratory tissue mask generation (Jupyter Notebook)
│   └── generate_masks.ipynb       # - tissue mask generation workflow
├── utils/                         # Helper functions for the workflows:
│   ├── canonical_markers.yaml     # - maps common marker synonyms to their canonical names
│   ├── io_utils.py                # - I/O helper functions
│   ├── mask_utils.py              # - tissue mask generation helper functions
│   ├── patch_utils.py             # - patch extraction helper functions
│   └── preprocessing.py           # - preprocessing helper functions
├── README.md
├── environment.yml
└── requirements.txt
```

## Installation

### Using Conda
```
git clone https://github.com/sophiamjiali/ws_imc_workflows.git
cd ws_imc_workflows
conda env create -f environment.yml
conda activate ws_imc_workflows_env
```

### Using pip
```
python -m pip install -r requirements.txt
```

## Configuration
All workflow parameters are stored in a YAML configuration file. A default configurations file is provided as `config/default.yaml` and a template provided as `config/template_config.yaml`. Configurations designed as required must be provided, where no default values are provided within the workflow. Those designated as not required will have default values provided if their corresponding optional process is enabled. 

### Using the configuration file
**1. Copy the template before editing**
```
cp config/template_config.yaml config/my_config.yaml
```
**2. Edit config/my_config.yaml to set appropriate paths and options.**

**3. Pass its path relative to the project root when running the script**
```
python -m src.generate_masks --config config/my_config.yaml
```

### Configuration sections
The YAML will contain the following sections with default values listed.

**1. General**
```
seed: 42                                          # Seed for reproducibility
input_folder: data/raw                            # Folder with raw/preprocessed images
panel_file: data/panels/ws_panel.csv              # Panel file for *all* provided images
id_mapping_file: data/id_mappings/id_mapping.csv  # Mapping file for internal WSI ID to image file paths
```
**2. Preprocessing**
```
preprocessing:
  toggles:
    apply_background_stain_removal: true          # Remove channels corresponding to background stains
    apply_hot_pixel_removal: true                 # Remove hot pixels per channel using a median filter
    apply_striping_removal: true                  # Remove striping artifacts per channel using a median filter
    apply_denoising: false                        # Denoises via a variance-stabilizing transform with arcsinh mapping
    apply_background_subtraction: false           # Subtracts per-channel background using a percentile estimate
    apply_winsorization: true                     # Clips the top and bottom quantile values per channel
    apply_min_max_scaling: true                   # Scales pixel intensities to [0, 1] per channel via min-max scaling

  background_stains:                              # Background markers to remove from the image(s) (can be left empty)
    - pnad                                        # Note: markers should be provided as their canonical marker name (^)

  hot_pixel:                                      # Parameters for hot pixel removal:
    window_size: 3                                # - Window size for the median filter
    z_score_threshold: 6                          # - Threshold for identifying candidate hot pixels by their z-score

  striping:                                       # Parameters for striping artifact removal:    
    direction: column                             # - Direction of striping artifacts (column/row) (manual inspection)
    size: 5                                       # - Window size for the median filter

  denoising:                                      # Parameters for denoising:
    cofactor: 10                                  # - Cofactor for arcsinh mapping per channel, per pixel

  background_subtraction:                         # Parameters for removing per-channel background:
    percentile: 2.0                               # - Percentile to estimate background intensity

  winsorization:                                  # Parameters for winsorizing per-channel intensity 
    limits: [0.005, 0.005]                        # - Top and bottom quantiles for clipping
```
**3. Tissue Mask**
```
tissue_mask:
  mask_folder: results/tissue_masks/masks         # Folder for output masks
  metadata_folder: results/tissue_masks/metadata  # Folder for mask metadata
  qc_folder: results/tissue_masks/qc              # Folder for mask overlay QC plots

  mask_generation_markers:                        # Markers to define the tissue mask upon (can be left empty)
    - pan_cytokeratin                             # Note: if no markers are provided, the tissue mask will be 
    - vimentin                                    #       generated on all available channels (after 
    - col1a2                                      #       preprocessing if enabled)
    - fibronectin
    - tenacin
    - cd45
    - acta2
    - fap
    - cd31

  min_tissue_threshold: 0                         # Minimum tissue threshold for dynamic threshold determination
  small_object_threshold: 50                      # Maximum size in pixels of an object to remove
  small_hole_threshold: 50                        # Maximum size in pixels of a hole to fill

  rgb_markers:                                    # Markers for RGB colour channels (in order) for QC plot visualization
    - pan_cytokeratin   # Red                     # Note: markers should be provided as their canonical marker name (see
    - cd45              # Green                   #       utils/canonical_markers.yaml) for reference
    - col1a2            # Blue
```
**4. Patch Extraction**
```
patch_extraction:
  patch_folder: results/patch_extraction/patches.zarr  # .zarr folder for patches, masks, and metadata
  mask_folder: results/tissue_masks/masks              # Folder with WSI masks (matched by WSI ID)
  patch_size: [128, 128]                               # Patch sizes for extraction
  stride: 0.25                                         # Stride for extraction
  min_tissue_coverage: 0.50                            # Minimum tissue coverage of each patch
```

## Usage

### 1. Running the workflow

The workflow is run from the project root using the following Python module syntax:

**A. Tissue Mask Generation**

This workflow will generate tissue masks for each IMC image provided using GMM/OTsu thresholding with optional preprocessing and cleaning.
```
python -m src.generate_masks \
  --config config/my_config.yaml \
  --threshold-method gmm \
  --remove-small-objects True \
  --fill-small-holes True \
  --save-metadata True \
  --save-qc True \
  --preprocess-images True
```
Output will include:
- binary tissue masks (`.tiff`)
- optional quality control visualizations (`.png`)
- optional mask generation metadata (`.csv`)

**B. Patch Extraction**

This workflow will extract patches for each IMC image provided with optional preprocessing.
```
python -m src.extract_patches --config config/my_config.yaml --preprocess-images True
```
Output will include:
- patches grouped by WSI ID (`.zarr`)
- patch tissue masks grouped by WSI ID (`.zarr`)
- manifest mapping patch ID, WSI ID, and metadata (`.csv`)
- patch extraction statistics (`.csv`)
- cohort extraction statistics (`.json`)

### 2. Input configuration

All inputs and options are specified in a YAML configuration file. See `config/template_config.yaml` and `config/default.yaml` for examples. A comprehensive breakdown of all required configurations is provided in both configuration files.

**Both Workflows**

You must at minimum specify:
- `seed`: seed for reproducibility
- `input_folder`: folder containing IMC `.mcd` or `.tiff` files
- `panel_file`: panel mapping markers to metal tags as a `.csv`
- `id_mapping_file`: mapping between internal WSI ID with image file paths

**A. Tissue Mask Generation**

You must at minimum specify:
- `tissue_mask.mask_folder`: folder for binary tissue masks

**B. Patch Extraction**

You must at minimum specify:
- `patch_extraction.patch_folder`: folder to store patches, masks, and metadata as a `.zarr`
- `patch_extraction.mask_folder`: folder containing binary tissue masks for *all* images provided
- `patch_extraction.patch_size`: list representing the height and width of the patches
- `patch_extraction.stride`: stride/overlap when extracting patches
- `patch_extraction.min_tissue_coverage`: minimum tissue coverage per patch

Non-required configurations corresponding to optional functionality have default values provided within the workflow if not provided.

### 3. Example folder structure

Note that the workflow can be performed on both raw and preprocessed IMC images. If the input is raw data and should be preprocessed prior to tissue mask generation, enable `--preprocess-images` to `True` when running `generate_masks.py` and set corresponding configurations in `config/my_config.yaml`.
```
ws_imc_workflows/
├── config/
│   ├── template_config.yaml
│   ├── default.yaml
│   └── my_config.yaml
├── data/
│   ├── raw/
│   │   ├── sample_01.mcd
│   │   ├── sample_02.mcd
│   │   └── ...
│   ├── preprocessed/
│   │   ├── sample_01.tiff
│   │   ├── sample_02.tiff
│   │   └── ...
│   └── panels/
│       ├── panel_raw_data.csv
│       ├── panel_processed_data.csv
│       └── ...
└── results/
    ├── tissue_masks/
    │   ├── masks/
    │   │   ├── sample_01_mask.tiff
    │   │   ├── sample_01_mask.tiff
    │   │   ├── sample_02_mask.tiff
    │   │   └── ...
    │   ├── metadata/
    │   │   ├── sample_01_mask_metadata.csv
    │   │   ├── sample_02_mask_metadata.csv
    │   │   └── ...
    │   └── qc/
    │       ├── sample_01_mask_qc.png
    │       ├── sample_02_mask_qc.png
    │       └── ...
    └── patch_extraction/
```

### 4. Output Files
After running the workflow, the following outputs will be generated:

### Tissue Mask Generation
Output            | File Type | Description                                     | Naming Convention                    
----------------- | --------- | ----------------------------------------------- | ------------------------------------ 
Image Masks       | `.tiff`   | Binary tissue masks per image                   | `{wsi_id}_mask.tiff`                 
Metadata          | `.csv`    | Threshold values, coverage %, method used, etc. | `{wsi_id}_mask_metadata.csv`         
QC Plots          | `.png`    | Composite and mask overlay visualization        | `{wsi_id}_mask_qc.png`               

### Patch Extraction
Output            | File Type | Description                                     | Naming Convention                    
----------------- | --------- | ----------------------------------------------- | ------------------------------------ 
Patches           | `.zarr`   | Patches extracted from preprocessed WSI         | `{wsi_id}_x{}_y{}_patch_{patch_idx`} 
Patch Masks       | `.json`   | Binary tissue masks per patch                   | `{wsi_id}_x{}_y{}_patch_{patch_idx`} 
Manifest          | `.csv`    | Patch ID, WSI ID, and patch metadata            | `manifest.csv`                       
Image Statistics  | `.csv`    | Number of attempted and valid patches per image | `wsi_statistics.csv`                 
Cohort Statistics | `.json`   | Cohort patch extraction statistics              | `cohort_statistics.json`
WSI ID Mapping    | `.csv`    | Mapping for internal WSI ID to image filename   | `id_mapping.csv`                     
