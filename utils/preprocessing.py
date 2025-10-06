"""
Script:         preprocessing.py
Purpose:        Functions for preprocessing images
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


def preprocess_image(image, config):
    # Returns the preprocessed xarray image per the configurations provided

    # Fetch thresholds from the configurations
    toggles_cfg = preproc_cfg.get('toggles', {})
    preproc_cfg = config.get('preprocessing', {})
    
    # 1. Background stain removal
    if toggles_cfg.get('apply_background_stain_removal', True):
        background_stains = preproc_cfg.get('background_stains', [])
        image = remove_background_stains(image, background_stains)

    # 2. Hot pixel removal
    if toggles_cfg.get('apply_hot_pixel_removal', True):
        hot_pixel_cfg = preproc_cfg.get('hot_pixel', {})
        image = remove_hot_pixels(image,
            window_size = hot_pixel_cfg.get('window_size', 0),
            z_score_threshold = hot_pixel_cfg.get('z_score_threshold', 0)
        )

    # 3. Striping artifact removal
    if toggles_cfg.get('apply_striping_removal', True):
        striping_cfg = preproc_cfg.get('striping', {})
        image = remove_striping_artifacts(image,
            direction = striping_cfg.get('direction', 'vertical'),
            size = striping_cfg.get('size', 0)
        )

    # 4. Shot-noise denoising
    if toggles_cfg.get('apply_denoising', True):
        denoise_cfg = preproc_cfg.get('denoising', {})
        image = denoise(image,
            cofactor = denoise_cfg.get('cofactor', 0)
        )

    # 5. Background subtraction
    if toggles_cfg.get('apply_background_subtraction', True):
        background_cfg = preproc_cfg.get('background_subtraction', {})
        image = subtract_background(image,
            percentile = background_cfg.get('percentile', 0)
        )

    # 6. Winsorization
    if toggles_cfg.get('apply_winsorization', True):
        winsorize_cfg = preproc_cfg.get('winsorization', {})
        image = winsorize(image, 
            limits = winsorize_cfg.get('limits', [0, 0])
        )

    # 7. Global Normalization
    if toggles_cfg.get('apply_min_max_scaling', True):
        image = scale(image)

    return image

def remove_background_stains(image, background_stains):
    # Removes the background stains, provided as metal tags, from the image
    return image.drop_sel(metal_tag = background_stains)

def remove_hot_pixels(image, window_size, z_score_threshold):
    # Removes hot pixels per channel using a median filter to estimate 
    # background, and removing only single hot pixel points; only non-isolated 
    # hot pixels are kept as candidates to avoid removing small clusters, which 
    # could be real biological features in 5 micron resolution
    
    # Define a hot pixel mask per channel
    clean_image = np.empty_like(image)

    # Define an 8-connected neighbourhood structure for each channel
    structure = generate_binary_structure(2, 1)

    # Detect and remove hot pixels independently per channel
    for c_idx in range(image.shape[0]):
        channel = image.isel(channel = c_idx)

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
        clean_channel[isolated_mask] = local_median[isolated_mask]

        clean_image[c_idx] = clean_channel

    return clean_image

def remove_striping_artifacts(image, direction, size):
    # Applies a directional median filter per channel along the given direction 
    # to reduce striping artifacts

    filtered_img = np.empty_like(image)

    if direction == "row":
        size_tuple = (1, int(size))
    elif direction == "column":
        size_tuple = (int(size), 1)
    else:
        raise ValueError("Direction must be 'row' or 'column'.")

    # Apply the median filter per channel
    for c in range(image.shape[0]):
        channel = np.ascontiguousarray(image[c], dtype=np.float32)
        filtered_img[c] = median_filter(
            channel, size = size_tuple, mode = "reflect"
        )

    return filtered_img

def denoise(image, cofactor):
    # Performs a variance-stabilizing transform (VST) on the IMC image using
    # arcsinh mapping per channel, per pixel

    # Convert the image to float32 for safe computation
    vst_img = image.astype(np.float32, copy = True)

    # Apply the arcsinh transform to each channel independently
    for c in range(image.shape[0]):
        vst_img[c] = np.arcsinh(image[c] / cofactor)

    return vst_img

def subtract_background(image, percentile):
    # Subtract a robust per-channel background offset using the specified
    # percentile as a conservative background estimate; intensities are clipped
    # at zero to avoid negatives

    bg_subtract_img = image.copy()

    for c in range(image.shape[0]):
        bg_level = np.percentile(image[c], percentile)
        bg_subtract_img[c] = np.clip(
            image[c] - bg_level, a_min = 0, a_max = None
        )

    return bg_subtract_img

def winsorize(image, limits):
    # Clips the top and bottom quantile values of the image
    
    wins_img = np.empty_like(image, dtype = np.float32)
    for c in range(image.shape[0]):

        ## Flatten the channel, winsorize, then reshape
        flat = image[c].flatten()
        win_flat = scipy_winsorize(flat, limits = limits)
        wins_img[c] = win_flat.reshape(image.shape[1:])

    return wins_img

def scale(image):
    # Scales the pixel intensities to [0, 1] per channel via min-max scaling
    
    scaled = image.astype(np.float32, copy = True)
    scaler = MinMaxScaler(feature_range = (0, 1))
    
    for c in range(image.shape[0]):
        flat = image[c].reshape(-1, 1)
        scaled_flat = scaler.fit_transform(flat)
        scaled[c] = scaled_flat.reshape(image.shape[1:])

    return scaled