from pathlib import Path
import yaml
import pandas as pd
from readimc import MCDFile
import tifffile as tf
import xarray as xr
import os
import re
import numpy as np


def load_input_paths(config):
    # Loads and returns a list of .mcd file paths

    # Verify the input folder exists
    input_folder = Path(config.get('input_folder', None))

    if input_folder is None:
        raise ValueError("Config must provide 'input_folder'")
    elif not input_folder.exists():
        raise FileNotFoundError(f"Input folder {input_folder} does not exist.")

    # Verify files exist in the folder and their type
    mcd_files = list(input_folder.glob("*.mcd"))
    tiff_files = list(input_folder.glob("*.tiff"))

    if tiff_files:
        return (mcd_files, "TIFF")
    elif mcd_files:
        return (mcd_files, "MCD")
    else:
        raise RuntimeError(f"No MCD or TIFF files found in {input_folder}")


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


def load_panel(mcd_file, file_type, config):
    # Loads and returns an IMC stain panel as a dataframe, filtered for markers 
    # present in the provided data and background stains if specified; assumes 
    # markers are uniformly present in all input imgs 

    # Verify the panel exists
    panel_path = config.get('panel_file', None)
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

    # Create a new column for canonical marker names
    panel['canonical_marker'] = canonicalize_markers(panel['marker'].tolist())

    # Sanitize the metal tags
    panel['canonical_metal_tag'] = (canonicalize_metal_tags(panel['metal_tag']
                                                            .tolist()))
    
    # Remove any metal tags not present in the provided input data
    img = load_image(mcd_file, file_type, panel)
    present_tags = img.coords['metal_tag'].values.tolist()
    panel = panel[panel['canonical_metal_tag'].isin(present_tags)]
    
    # Filter out background stains (markers) if provided
    background_stains = config.get('background_stains', [])
    if not background_stains:
        panel = panel[~panel['marker'].isin(background_stains)]

    # Collapse any duplicate metal tags to the first channel
    panel = panel.drop_duplicates(subset = 'metal_tag', keep = 'first')

    return panel


def canonicalize_markers(markers):
    # Maps the given list of markers to its canonical name

    # Load the canonical markers mapping file
    mapping_path = Path('utils/canonical_markers.yaml')
    if not mapping_path.exists():
        raise FileNotFoundError(
            f"YAML Canonical Marker mapping file {mapping_path} does not exist"
        )
    try:
        with open(mapping_path, 'r') as f:
            mapping = yaml.safe_load(f)
    except mapping_path.YAMLError as e:
        raise RuntimeError(
            f"Failed to parse YAML canonical marker mapping file " +
            f"{mapping_path}: {e}"
        )
    
    # Build the reverse lookup: map display + synonyms + key → key
    value_to_key = {}
    for key, info in mapping.items():
        # Map the key itself
        value_to_key[key] = key
        # Map the display name
        display_name = info.get("display")
        if display_name is not None:
            value_to_key[display_name] = key
        # Map synonyms
        for syn in info.get("synonyms", []):
            value_to_key[syn] = key

    # Map the list of markers to their canonical names
    canonical_markers = [value_to_key.get(m, m) for m in markers]
    return canonical_markers


def canonicalize_metal_tags(metal_tags):
    # Sanitizes and returns the given metal tag 

    canonicalized = []
    for tag in metal_tags:
        if not isinstance(tag, str):
            tag = str(tag)
        tag = tag.replace(' ', '')
        letters = ''.join(re.findall(r'[A-Za-z]', tag))
        numbers = ''.join(re.findall(r'\d', tag))
        canonicalized.append(f"{letters}{numbers}")
    return canonicalized

           
def load_image(img_path, file_type, panel):
    # Loads and returns an image of the given file type as an xarray object

    if file_type == "MCD":

        # Verify the image can be parsed
        try:
            with MCDFile(img_path) as f:
                acquisition = f.slides[0].acquisitions[0]
                img = f.read_acquisition(acquisition)
                metal_tags = acquisition.channel_names
        except Exception as e:
            print(f"Error processing file: {img_path} - {e}")

    else:

        # Assuming the panel is one-to-one and index matched with the image
        img = tf.imread(img_path)
        metal_tags = panel['metal_tag'].tolist()

    # Return the img, markers, and metal tag labels as an xarray
    metal_tags = canonicalize_metal_tags(metal_tags)
    img = xr.DataArray(
        img, 
        dims = ("channel", "y", "x"),
        coords = {'metal_tag': ('channel', metal_tags)}
    )  

    return img

def to_xarray(data, template):
    # Returns an xarray using the provided data with the template's metadata
    return xr.DataArray(
        data, 
        dims = template.dims, 
        coords = template.coords, 
        attrs = template.attrs
    )


def load_mask_generation_panel(config, panel):
    # Loads and returns dataframe of metal tag - marker mappings to be used to 
    # generate the tissue mask upon; if a subset was not provided, the full 
    # panel is returned

    # Fetch the tissue mask generation markers if provided
    mask_panel = panel.copy()
    mask_markers = (config.get('tissue_mask', {})
                          .get('mask_generation_markers', []))

    # Subset for the provided mask markers, else use all markers present
    if mask_markers:
        mask_panel = mask_panel[
            mask_panel['canonical_marker'].isin(mask_markers)
        ]
    
    return mask_panel


def save_tissue_mask(tissue_mask, img_name, config):
    # Saves the provided tissue mask using the name of the img with a 
    # '_mask.tiff' suffix to the output folder specified in the configurations

    # If the output folder was not initialized, create it using its path
    mask_folder = Path(config.get('tissue_mask', {})
                           .get('mask_folder', None))

    if mask_folder is None:
        raise ValueError("Config must provide 'tissue_mask.mask_folder'")
    mask_folder.mkdir(parents = True, exist_ok = True)
        
    out_path = os.path.join(mask_folder, img_name + "_mask.tiff")
    tf.imwrite(out_path, tissue_mask.astype(np.uint8))


def save_mask_metadata(metadata, img_name, config):
    # Saves the provided metadata as a CSV using the name of the img with 
    # a '_mask_metadata.csv' suffix to the output folder specified in the 
    # configurations

    # If a metadata folder was not provided, create one
    metadata_folder = Path(config.get('tissue_mask', {})
                                 .get('metadata_folder', None))
    
    if metadata_folder is None:
        raise ValueError("Config must provide 'tissue_mask.metadata_folder'")
    metadata_folder.mkdir(parents = True, exist_ok = True)
        
    # Save the metadata as a CSV in the folder
    out_path = os.path.join(metadata_folder, img_name + "_metadata.csv")
    metadata = pd.DataFrame([metadata])
    metadata.to_csv(out_path, index = False)
        

def save_mask_qc(qc_plot, img_name, config):
    # Saves the provided quality control plot as a PNG using the name of the
    # img with a '_qc.png' suffix to the output folder specified in the
    # configurations

    # If a quality control folder was not provided, create one
    qc_folder = Path(config.get('tissue_mask', {}).get('qc_folder', None))

    if qc_folder is None:
        raise ValueError("Config must provide 'tissue_mask.qc_folder'")
    qc_folder.mkdir(parents = True, exist_ok = True)
        
    # Save the QC plot as a PNG in the folder
    out_path = os.path.join(qc_folder, img_name + "_qc.png")
    qc_plot.savefig(out_path, dpi = 300, bbox_inches = 'tight')


def generate_wsi_id_mapping(img_files, config):
    # Generates, returns, and saves a mapping from WSI file paths to internal 
    # WSI IDs

    # If a mapping ID folder was not provided, create one
    id_mapping_file = Path(config.get('id_mapping_file', None))
    if id_mapping_file is None:
        raise ValueError("Config must provide 'id_mapping_file'")
    id_mapping_file.parent.mkdir(parents = True, exist_ok = True)

    # Sort the files before mapping them to IDs for reproducibility
    img_files = sorted(img_files)
    file_mapper = [
        {'wsi_id': f"wsi_{i}", 'file_path': str(f)}
        for i, f in enumerate(img_files)
    ]
    
    # Save the mapping file as a CSV
    pd.DataFrame(file_mapper).to_csv(id_mapping_file, index = False)
    return file_mapper