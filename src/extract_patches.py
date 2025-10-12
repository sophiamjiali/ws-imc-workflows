#!/usr/bin/env python3
"""
Script:          extract_patches.py
Purpose:         Extracts patches, metadata, and statistics
Author:          Sophia Li
Affiliation:     Campbell Lab
Date:            010-03-2025
"""

import shutil
import argparse
import os
import pandas as pd
import zarr
import json
from pathlib import Path

from utils.patch_utils import extract_patches
from utils.io_utils import (
    load_input_paths,
    load_image,
    load_panel,
    load_config,
    generate_wsi_id_mapping
)

def parse_arguments():
    # Parses and returns command-line arguments

    parser = argparse.ArgumentParser(
        description = "Extraction patches from WS-IMC images."
    )
    parser.add_argument("--config", type = str, required = True,
                        help = "Path to configuration file")
    parser.add_argument('--preprocess-images', type = bool, default = True,
                        help = "Preprocess the images in the input folder")
    
    return parser.parse_args()


def main():

    args = parse_arguments()
    config = load_config(args.config)
    img_files, file_type = load_input_paths(config)
    id_mapper = generate_wsi_id_mapping(img_files, config)
    panel = load_panel(img_files[0], file_type, config)

    # Fetch configurations
    patch_folder = Path(config.get('patch_extraction', {})
                              .get('patch_folder', None))
    if patch_folder is None:
        raise ValueError("Config must provide 'patch_extraction.patch_folder'")
    
    # Verify that tissue masks have been generated
    mask_folder = Path(config.get('patch_extraction', {})
                             .get('mask_folder', None))
    if mask_folder is None:
        raise ValueError(f"Mask folder {mask_folder} does not exist. " +
                         "If tissue masks have not been generated, " + 
                         "perform the tissue mask generation " +
                         "workflow+ prior to extracting patches.")

    # Fetch the number of channels per image (assumed the same for all)
    test_img = load_image(img_files[0], file_type, panel)
    C = test_img.shape[0]
    H, W = config.get('patch_extraction', {}).get('patch_size', [0, 0])

    # Initialize the OME.zarr folder to record the patches and metadata
    zarr_root = zarr.open(patch_folder, mode = 'w')
    patch_group = zarr_root.create_group('patches')
    mask_group = zarr_root.create_group('masks')

    # Record global metadata
    patch_idx = 0
    manifest, patch_statistics = [], []
    cohort_stats = {'total_attempted_patches': 0, 'total_valid_patches': 0}

    try:
        for path in img_files:

            # Map the filepath to an internal WSI ID and create subgroups
            wsi_id = id_mapper[str(path)]
            img_patch_group = patch_group.create_group(wsi_id)
            img_mask_group = mask_group.create_group(wsi_id)

            # Enter per patch into the .zarr to avoid memory overhead
            stats = {'total_attempted': 0, 'total_valid': 0}

            for patch, p_mask, meta in extract_patches(path, file_type, wsi_id,
                                                    panel, stats, config, 
                                                    args.preprocess_images):
                
                # Record the patch into the .zarr and build up the patch manifest
                name = f"{wsi_id}_y{meta['y']}_x{meta['x']}_patch_{patch_idx}"
                img_patch_group[name] = patch
                img_mask_group[name] = p_mask
                manifest.append({'patch_idx': patch_idx, **meta})

                patch_idx += 1

            # Update global metadata after each image is processed
            patch_statistics.append({
                'wsi_id': wsi_id, 
                'attempted_patches': stats['total_attempted'],
                'valid_patches': stats['total_valid']
            })
            cohort_stats['total_attempted_patches'] += stats['total_attempted']
            cohort_stats['total_valid_patches'] += stats['total_valid']

        # Write the metadata as CSVs and a JSON
        manifest_path = Path(os.path.join(
            patch_folder, 'metadata', 'manifest.csv'
        ))
        manifest_path.parent.mkdir(parents = True, exist_ok = True)
        pd.DataFrame(manifest).to_csv(manifest_path, index = False)

        stats_path = Path(os.path.join(
            patch_folder, 'metadata', 'wsi_statistics.csv'
        ))
        pd.DataFrame(patch_statistics).to_csv(stats_path, index = False)

        cohort_stats_path = Path(os.path.join(
            patch_folder, 'metadata', 'cohort_statistics.json'
        ))
        with open(cohort_stats_path, 'w') as f:
            json.dump(cohort_stats, f, indent = 2)

    except:
        shutil.rmtree(patch_folder)

if __name__ == '__main__':
    main()