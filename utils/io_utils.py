from pathlib import Path
import yaml
import pandas as pd
from readimc import MCDFile
import xarray as xr
import os

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
    mask_markers = config.get('mask_generation_markers', [])

    # Subset for the provided mask markers, else use all markers present
    if mask_markers:
        mask_panel = mask_panel[mask_panel['marker'].isin(mask_markers)]
    
    return mask_panel

def save_tissue_mask(tissue_mask, image_name, config):
    # Saves the provided tissue mask using the name of the image with a 
    # '_mask.tiff' suffix to the output folder specified in the configurations
    return

def save_mask_metadata(metadata, image_name, config):
    # Saves the provided metadata as a CSV using the name of the image with 
    # a '_mask_metadata.csv' suffix to the output folder specified in the 
    # configurations

    # If a metadata folder was not provided, create one
    metadata_folder = config.get('paths', {}).get('metadata_folder', None)

    if not metadata_folder:
        metadata_folder = os.path.join(
            os.path.join(config['paths']['output_folder']), 'metadata'
        ).mkdir(parents = True, exist_ok = True)

        print("No metadata folder specified in the configuration. Using" +
              f"default {metadata_folder}.")
        
    # Save the metadata as a CSV in the folder
    out_path = os.path.join(metadata_folder, image_name + "_metadata.csv")
    metadata.to_csv(out_path, index = False)
        

def save_mask_qc(qc_plot, image_name, config):
    # Saves the provided quality control plot as a PNG using the name of the
    # image with a '_qc.png' suffix to the output folder specified in the
    # configurations

    # If a quality control folder was not provided, create one
    qc_folder = config.get('paths', {}).get('qc_folder', None)

    if not qc_folder:
        qc_folder = os.path.join(
            os.path.join(config['paths']['output_folder']), 'qc'
        ).mkdir(parents = True, exist_ok = True)

        print("No quality control (qc) folder specified in teh configuration." +
              f" Using default {qc_folder}.")
        
    # Save the QC plot as a PNG in the folder
    out_path = os.path.join(qc_folder, image_name + "_qc.png")
    image_name ...