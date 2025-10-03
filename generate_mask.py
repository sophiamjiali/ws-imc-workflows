#!/usr/bin/env python3

import re
import numpy as np
import pandas as pd
import argparse
import yaml
import os
from pathlib import Path
import sys
import xarray as xr
from readimc import MCDFile

from sklearn.mixture import GaussianMixture
from skimage.filters import threshold_otsu
from skimage import morphology, measure, exposure
from skimage.measure import find_contours
from skimage.morphology import remove_small_objects, remove_small_holes 

def load_input_paths(input_path):
    # Loads and returns a list of .mcd file paths

    # Verify the input folder exists
    input_folder = Path(input_path)
    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder {input_folder} does not exist.")
    
    # Verify .mcd files exist in the folder
    mcd_files = list(input_folder.glob("*.mcd"))
    if not mcd_files:
        raise RuntimeError(f"No .mcd files found in {input_folder}")
    
    return mcd_files

def load_config(config_path):
    # Loads and returns a YAML configuration file

    # Verify the configuration file exists
    config = Path(config_path)
    if not config.exists():
        raise FileNotFoundError(
            f"YAML configuration file {config} does not exist"
        )

    # Verify the configuration file can be parsed
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except config.YAMLError as e:
        raise RuntimeError(
            f"Failed to parse YAML configuration file {config}: {e}"
        )
    
    return config

def load_panel(panel_path, mcd_file, config):
    # Loads and returns an IMC stain panel as a dataframe, filtered for markers 
    # present in the provided data and background stains if specified; assumes 
    # markers are uniformly present in all input images 

    # Verify the panel exists
    panel_file = Path(panel_path)
    if not panel_file.exists():
        raise FileNotFoundError(f"Panel {panel_file} does not exist")

    # Verify the configuration file can be parsed
    if panel_file.suffix.lower() == '.csv':
        try:
            panel = pd.read_csv(panel_file)
        except Exception as e:
            raise RuntimeError(f"Failed to read CSV panel {panel_file}: {e}")
    elif panel_file.suffix.lower() in ['.xls', '.xlsx']:
        try:
            panel = pd.read_excel(panel_file)
        except Exception as e:
            raise RuntimeError(f"Failed to read Excel panel {panel_file}: {e}")
    else:
        raise ValueError(f"Unsupported panel file format: {panel_file.suffix}")
    
    # Identify the metal tag and marker columns and standardize their names
    metal_col = next((c for c in ['Metal', 'MetalTag'] if c in panel.columns))
    marker_col = next((c for c in ['Target', 'Marker'] if c in panel.columns))
    panel = panel[[metal_col, marker_col]]
    panel.columns = ['metal_tag', 'marker']
    
    # Standardize value types to strings
    panel['metal_tag'] = panel['metal_tag'].astype(str)

    # Filter out background stains (markers) if provided
    background_stains = config.get('background_stains, []')
    if not background_stains:
        panel = panel[~panel['marker'].isin(background_stains)]

    # Remove any metal tags not present in the provided input data
    image = load_image(mcd_file)
    present_tags = list(image.metal_tag.values)
    panel = panel[panel['metal_tag'].isin(present_tags)]

    return panel
           
def load_image(image_path):
    # Loads and returns an .mcd image as an xarray object

    # Verify the .mcd image can be parsed
    try:
        with MCDFile(image_path) as f:
            acquisition = f.slides[0].acquisitions[0]
            image = f.read_acquisition(acquisition)
            metal_tags = acquisition.channel_names
    except Exception as e:
        print(f"Error processing file: {image_path} - {e}")

    # Return the image and metal tag labels as an xarray
    image = xr.DataArray(
        image, 
        dims = ("channel", "y", "x"),
        coords = {'metal_tag': metal_tags}
    )

    return image

def load_mask_generation_panel(config, panel):
    # Loads and returns dataframe of metal tag - marker mappings to be used to 
    # generate the tissue mask upon; if a subset was not provided, the full 
    # panel is returned

    # Fetch the tissue mask generation markers if provided
    mask_panel = panel.copy()

    mask_markers = config.get('mask_generation_markers, []')
    if mask_markers:
        mask_panel = mask_panel[mask_panel['marker'].isin(mask_markers)]
    
    return mask_panel


def generate_mask():
    return

def parse_arguments():
    # Parses and returns command-line arguments

    parser = argparse.ArgumentParser(
        description = "Generate tissue masks from WS-IMC images."
    )
    parser.add_argument("--input", type = str, required = True, 
                        help = "Path to input folder")
    parser.add_argument("--output", type = str, required = True,
                        help = "Path to output folder")
    parser.add_argument("--config", type = str, default = 'config.yaml',
                        help = "Path to configuration file")
    parser.add_argument("--panel", type = str, required = True,
                        help = "Path to panel")
    parser.add_argument("--threshold-method", type = str, 
                        choices = ['otsu', 'gmm'], default = 'otsu'
                        help = "Thresholding method to determine a tissue threshold for the tissue mask: 'otsu', or 'gmm'")
    parser.add_argument("--remove-small-objects", type = bool, default = True,
                        help = "Remove small objects from the tissue mask")
    parser.add_argument("--fill-small-holes", type = bool, default = True,
                        help = "Fill small holes in the tissue mask")
    return parser.parse_args()


def main():

    # Parse command-line arguments
    args = parse_arguments()
    
    # Load configurations
    config = load_config(args.config)

    # Load input paths
    mcd_files = load_input_paths(args.input)

    # Load the panel
    panel = load_panel(args.panel, mcd_files[0], config)

    # Fetch mask generation markers if provided
    mask_panel = load_mask_generation_panel(config, panel)
    mask_metals = mask_panel['metal_tag'].tolist()

    # Generate a tissue mask for each input image
    for path in mcd_files:

        # Load the image and filter for markers used for mask generation
        image = load_image(path)
        composite = image.sel(metal_tag = mask_metals)

        # Dynamically generate a tissue threshold for the image
        if args.threshold_method == 'otsu':
            threshold = determine_otsu_tissue_threshold(composite, config)
        else: # GMM
            threshold = determine_gmm_tissue_threshold(composite, config)


def determine_otsu_tissue_threshold(composite, config):
    # Determines and returns a global intensity threshold to separate tissue 
    # from background using Otsu's method

    # Fetch thresholds from the configurations
    clip_percentile = config.get('clip_percentile')

    # Clip extreme pixels of high intensity per channel
    clipped_channels = np.array([
        np.minimum(
            composite.sel(metal_tag = metal_tag).values,
            np.percentile(
                composite.sel(metal_tag, metal_tag).values, clip_percentile
            )
        ) for metal_tag in composite.metal_tag.values
    ])
    return

def determine_gmm_tissue_threshold(composite, config):
    # Determines and returns a global intensity threshold to separate tissue 
    # from background using a 2-component GaussianMixture Model (GMM) approach
