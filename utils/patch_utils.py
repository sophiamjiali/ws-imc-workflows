"""
Script:         patch_utils.py
Purpose:        Functions for patch extraction and manifest creation
Author:         Sophia Li
Affiliation:    Campbell Lab, Lunenfeld-Tanenbaum Research Institute (LTRI),
                University of Toronto
Date:           October 11 2025
"""

from utils.io_utils import (
    load_image
)
import os
import tifffile as tf
import numpy as np
from utils.preprocessing import preprocess_image


def extract_patches(img_path, file_type, wsi_id, panel, 
                    stats, config, preprocess):
    # Extracts and returns patches from the specified image, screening for 
    # sufficient tissue coverage; patches that overlap the image boundary
    # are discarded

    mask_folder = config.get('patch_extraction', {}).get('mask_folder', None)
    patch_size = config.get('patch_extraction', {}).get('patch_size', [0, 0])
    min_coverage = (config.get('patch_extraction', {})
                          .get('min_tissue_coverage', 0))
    stride = config.get('patch_extraction', {}).get('stride', 1.0)
    stride_H, stride_W = int(patch_size[0]*stride), int(patch_size[1]*stride)

    # Extract all potential patches and extract their metadata
    img = load_image(img_path, file_type, panel)
    if preprocess:
        img = preprocess_image(img, config)
    H, W = img.shape[-2:]

    # Find the corresponding tissue mask and load it, based on the image's name
    mask_file = os.path.join(mask_folder, wsi_id + '_mask.tiff')
    if os.path.exists(mask_file):
        mask = tf.imread(mask_file)
    else:
        raise FileNotFoundError(f"Mask for image {img_path} does not exist.")
    
    total_attempted, total_valid = 0, 0

    # Compute and QC patches using a sliding window approach
    for y in range(0, H - patch_size[0] + 1, stride_H):
        for x in range(0, W - patch_size[1] + 1, stride_W):
            total_attempted += 1
            
            # Extract the patch and its corresponding mask
            patch = np.asarray(img.isel(
                y = slice(y, y + patch_size[0]), 
                x = slice(x, x + patch_size[1])
            ))
            patch_mask = mask[y:y + patch_size[0], x:x +patch_size[1]]

            # Screen for sufificent tissue coverage
            if np.mean(patch_mask) < min_coverage:
                continue

            total_valid += 1 

            # Return the patch and its metadata
            metadata = {
                'wsi_id': wsi_id,
                'y': y,
                'x': x,
                'stride': stride,
                'channels': patch.shape[0],
                'height': patch.shape[1],
                'width': patch.shape[2],
                'coverage': patch_mask.mean()
            }

            yield patch, patch_mask, metadata

    # Update mutable stats to return
    stats['total_attempted'] += total_attempted
    stats['total_valid'] += total_valid
