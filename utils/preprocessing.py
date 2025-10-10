"""
Script:         preprocessing.py
Purpose:        Functions for preprocessing imgs
Author:         Sophia Li
Affiliation:    Campbell Lab, Lunenfeld-Tanenbaum Research Institute (LTRI),
                University of Toronto
Date:           October 5 2025
"""

import numpy as np
from scipy.stats.mstats import winsorize as scipy_winsorize
from sklearn.preprocessing import MinMaxScaler
from scipy.ndimage import (
    median_filter, 
    generate_binary_structure,
    binary_erosion
)
from utils.io_utils import to_xarray


def preprocess_image(img, config):
    # Returns the preprocessed xarray img per the configurations provided

    # Fetch thresholds from the configurations
    preproc_cfg = config.get('preprocessing', {})
    toggles_cfg = preproc_cfg.get('toggles', {})
    
    # 1. Background stain removal
    if toggles_cfg.get('apply_background_stain_removal', True):
        background_stains = preproc_cfg.get('background_stains', [])
        img = remove_background_stains(img, background_stains)

    # 2. Hot pixel removal
    if toggles_cfg.get('apply_hot_pixel_removal', True):
        hot_pixel_cfg = preproc_cfg.get('hot_pixel', {})
        processed_img = remove_hot_pixels(img,
            window_size = hot_pixel_cfg.get('window_size', 0),
            z_score_threshold = hot_pixel_cfg.get('z_score_threshold', 0)
        )
        img = to_xarray(processed_img, img)

    # 3. Striping artifact removal
    if toggles_cfg.get('apply_striping_removal', True):
        striping_cfg = preproc_cfg.get('striping', {})
        processed_img = remove_striping_artifacts(img,
            direction = striping_cfg.get('direction', 'vertical'),
            size = striping_cfg.get('size', 0)
        )

        img = to_xarray(processed_img, img)

    # 4. Shot-noise denoising
    if toggles_cfg.get('apply_denoising', True):
        denoise_cfg = preproc_cfg.get('denoising', {})
        processed_img = denoise(img,
            cofactor = denoise_cfg.get('cofactor', 0)
        )
        img = to_xarray(processed_img, img)

    # 5. Background subtraction
    if toggles_cfg.get('apply_background_subtraction', True):
        background_cfg = preproc_cfg.get('background_subtraction', {})
        processed_img = subtract_background(img,
            percentile = background_cfg.get('percentile', 0)
        )
        img = to_xarray(processed_img, img)

    # 6. Winsorization
    if toggles_cfg.get('apply_winsorization', True):
        winsorize_cfg = preproc_cfg.get('winsorization', {})
        processed_img = winsorize(img, 
            limits = winsorize_cfg.get('limits', [0, 0])
        )
        img = to_xarray(processed_img, img)

    # 7. Global Normalization
    if toggles_cfg.get('apply_min_max_scaling', True):
        processed_img = scale(img)
        img = to_xarray(processed_img, img)

    return img

def remove_background_stains(img, background_stains):
    # Removes the background stains, provided as metal tags, from the img
    return img.sel(channel = ~np.isin(img['metal_tag'], background_stains))

def remove_hot_pixels(img, window_size, z_score_threshold):
    # Removes hot pixels per channel using a median filter to estimate 
    # background, and removing only single hot pixel points; only non-isolated 
    # hot pixels are kept as candidates to avoid removing small clusters, which 
    # could be real biological features in 5 micron resolution
    
    # Define a hot pixel mask per channel
    clean_img = np.empty_like(img)

    # Define an 8-connected neighbourhood structure for each channel
    structure = generate_binary_structure(2, 1)

    # Detect and remove hot pixels independently per channel
    for c_idx in range(img.shape[0]):
        channel = img.isel(channel = c_idx)

        # Compute the local median and local median absolute deviation (MAD)
        local_median = median_filter(channel, window_size)
        local_mad = median_filter(np.abs(channel - local_median), window_size)

        # Pad regions with constant intensity to prevent zero division
        local_mad[local_mad == 0] = 1e-6

        # Calculate each pixel's z-score: how many MADs away from local median
        z_score = np.divide(
            (channel - local_median), 
            local_mad,
            out = np.zeros_like(channel, dtype = np.float32),
            where = local_mad != 0
        )

        # Identify candidate hot pixels based on the z-score threshold
        hot_pixel_mask = np.abs(z_score) > z_score_threshold

        # Restrict to only isolated pixels for hot pixel candidates
        isolated_mask = hot_pixel_mask & ~binary_erosion(
            hot_pixel_mask, structure
        )

        # Replace the isolated hot pixels with the local median value
        clean_channel = channel.copy()
        clean_channel = (clean_channel.where(
            ~isolated_mask, other = local_median
        ))

        clean_img[c_idx] = clean_channel

    return clean_img

def remove_striping_artifacts(img, direction, size):
    # Applies a directional median filter per channel along the given direction 
    # to reduce striping artifacts

    filtered_img = np.empty_like(img)

    if direction == "row":
        size_tuple = (1, int(size))
    elif direction == "column":
        size_tuple = (int(size), 1)
    else:
        raise ValueError("Direction must be 'row' or 'column'.")

    # Apply the median filter per channel
    for c in range(img.shape[0]):
        channel = np.ascontiguousarray(img[c], dtype=np.float32)
        filtered_img[c] = median_filter(
            channel, size = size_tuple, mode = "reflect"
        )

    return filtered_img

def denoise(img, cofactor):
    # Performs a variance-stabilizing transform (VST) on the IMC img using
    # arcsinh mapping per channel, per pixel

    # Convert the img to float32 for safe computation
    vst_img = img.astype(np.float32, copy = True)

    # Apply the arcsinh transform to each channel independently
    for c in range(img.shape[0]):
        vst_img[c] = np.arcsinh(img[c] / cofactor)

    return vst_img

def subtract_background(img, percentile):
    # Subtract a robust per-channel background offset using the specified
    # percentile as a conservative background estimate; intensities are clipped
    # at zero to avoid negatives

    bg_subtract_img = img.copy()

    for c in range(img.shape[0]):
        bg_level = np.percentile(img[c], percentile)
        bg_subtract_img[c] = np.clip(
            img[c] - bg_level, a_min = 0, a_max = None
        )

    return bg_subtract_img.values

def winsorize(img, limits):
    # Clips the top and bottom quantile values of the img
    
    wins_img = np.empty_like(img.values, dtype = np.float32)
    for c in range(img.shape[0]):

        ## Flatten the channel, winsorize, then reshape
        flat = img[c].values.flatten()
        win_flat = scipy_winsorize(flat, limits = limits)
        wins_img[c] = win_flat.reshape(img.shape[1:])

    return wins_img

def scale(img):
    # Scales the pixel intensities to [0, 1] per channel via min-max scaling
    
    scaled = np.empty_like(img.values, dtype = np.float32)
    scaler = MinMaxScaler(feature_range = (0, 1))
    
    for c in range(img.shape[0]):
        flat = img[c].values.reshape(-1, 1)
        scaled_flat = scaler.fit_transform(flat)
        scaled[c] = scaled_flat.reshape(img.shape[1:])

    return scaled