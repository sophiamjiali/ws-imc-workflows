"""
Script:         mask_utils.py
Purpose:        Functions for mask computation and processing
Author:         Sophia Li
Affiliation:    Campbell Lab, Lunenfeld-Tanenbaum Research Institute (LTRI),
                University of Toronto
Date:           October 5 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from skimage.filters import threshold_otsu
from sklearn.exceptions import ConvergenceWarning
from skimage.morphology import remove_small_objects, remove_small_holes 


def determine_otsu_tissue_threshold(composite, config):
    # Determines and returns a global intensity threshold to separate tissue 
    # from background using Otsu's method

    # Fetch configurations
    min_tissue_threshold = (config
                            .get('preprocessing', {})
                            .get('min_tissue_threshold', 0))
    
    # Generate the threshold using Otsu's method and enforce a minimum
    threshold = threshold_otsu(composite)
    threshold = max(threshold, min_tissue_threshold)
    metadata = {
        'method': 'otsu',
        'threshold': threshold,
        'min_threshold': min_tissue_threshold,
    }

    return threshold, metadata


def determine_gmm_tissue_threshold(composite, config):
    # Determines and returns a global intensity threshold to separate tissue 
    # from background using a 2-component GaussianMixture Model (GMM) approach

    # Fetch configurations
    min_tissue_threshold = (config
                            .get('preprocessing', {})
                            .get('min_tissue_threshold', 0))
    seed = config.get('seed', 42)

    # Aggregate across channels using the median and flatten to pixels
    composite = np.median(composite, axis = 0)
    pixels = composite.flatten()

    try:
        # Fit a 2-component GaussianMixture Model (GMM)
        gmm = GaussianMixture(n_components = 2, random_state = seed)
        gmm.fit(pixels)

        # Find the intersection between the two Gaussians, sorting w.r.t. means
        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_).flatten()
        weights = gmm.weights_.flatten()

        order = np.argsort(means)
        means, stds, weights = means[order], stds[order], weights[order]

        threshold = find_gaussian_intersection(
            means[0], stds[0], weights[0], means[1], stds[1], weights[1]
        )

        # Enforce a minimum tissue threshold safeguard upon the intersection
        threshold = max(threshold, min_tissue_threshold)
        metadata = {
            'method': 'gmm',
            'threshold': threshold,
            'min_threshold': min_tissue_threshold,
            'gmm_means': means,
            'gmm_stds': stds,
            'gmm_weights': weights
        }
    
    except(ConvergenceWarning, ValueError) as e:
        print(f"GMM failed: {e}. Using Otsu fallback.")
        threshold, metadata = determine_otsu_tissue_threshold(composite, config)

    return threshold, metadata


def find_gaussian_intersection(m0, s0, w0, m1, s1, w1):
    # Returns the intersection between two Gaussians, provided their means,
    # standard deviations, and weights

    # Handle identical variances to avoid dividing by zero
    if np.isclose(s0, s1):
        return (m0 + m1) / 2
    
    # Coefficients for the quadratic equation a*x^2 + b*x + c = 0
    a = s0**2 - s1**2
    b = 2 * (m0 * s1**2 - m1 * s0**2)
    c = (m1**2 * s0**2 - m0**2 * s1**2
         + 2 * s0**2 * s1**2 * np.log((s1 * w0) / (s0 * w1)))

    # Discriminant
    disc = b**2 - 4 * a * c

    if disc < 0:
        # No real intersection, fallback to midpoint
        return (m0 + m1) / 2

    # Compute both solutions
    x1 = (-b + np.sqrt(disc)) / (2 * a)
    x2 = (-b - np.sqrt(disc)) / (2 * a)

    # Return the value between the two means
    return x1 if m0 < x1 < m1 else x2


def generate_tissue_mask(composite, threshold, remove_objects, 
                         fill_holes, config):
    # Generates and returns a binary tissue mask upon the composite, removing
    # artifacts based on configuration specifications

    # Fetch configurations
    mask_cfg = config.get('tissue_mask')
    object_threshold = mask_cfg.get('small_object_threshold')
    hole_threshold = mask_cfg.get('small_hole_threshold')

    # Generate the (H, W) mask and record the original area
    mask = (composite > threshold).any(dim = 'marker').values
    orig_area_px = np.sum(mask)

    original_mask = mask.copy()

    # Remove small objects if toggled, measuring them
    if remove_objects:
        mask = remove_small_objects(mask, object_threshold)
        removed_objects = original_mask & ~mask
        total_removed_area = np.sum(removed_objects)

    # Fill small holes if toggled, measuring them
    if fill_holes:
        mask = remove_small_holes(mask, hole_threshold)
        filled_holes = mask & ~original_mask
        total_filled_area = np.sum(filled_holes)

    # Record mask metadata
    final_area_px = np.sum(mask)
    img_area_px = composite.shape[1] * composite.shape[2]

    metadata = {
        'raw_mask_area_px': int(orig_area_px),
        'clean_mask_area_px': int(final_area_px),
        'full_img_area_px': int(img_area_px),
        "removed_objects_area_px": int(total_removed_area),
        "filled_holes_area_px": int(total_filled_area),
        "mask_coverage_percent": int(final_area_px / img_area_px),
        "cleaning_area_change_px": int(final_area_px - orig_area_px),
        "cleaning_area_change_percent": ((final_area_px - orig_area_px) 
                                / orig_area_px * 100 
                                if orig_area_px > 0 else 0)
    }

    return mask, metadata


def generate_mask_qc_plot(mask, composite, image_name, 
                          threshold_metadata, mask_metadata, config):
    # Generates and returns a quality control (QC) plot overlaying the tissue
    # mask upon the composite image

    # Filter the composite for the RGB markers specified in the configurations
    rgb_markers = config.get('tissue_mask', {}).get('rgb_markers')
    composite = composite.sel(marker = rgb_markers)
    composite.transpose('y', 'x', 'marker')

    # Plot and return the figure
    fig, ax = plt.subplots(figsize = (8, 8))
    ax.imshow(composite, cmap = 'gray')
    ax.imshow(mask, cmap = 'Reds', alpha = 0.3)
    
    plt.set_title(f"Tissue Mask QC - {image_name}")

    # Build a caption based on the threshold method metadata
    ax.text(0.02, 0.98, 
                f"Method: GMM\nThreshold: {threshold_metadata['threshold']}\nCoverage: {mask_metadata['mask_coverage_percent']}",
                fontsize = 10, color = 'white', transform = plt.gca().transAxes, verticalalignment = 'top', bbox = dict(facecolor = 'black', alpha = 0.5)
    )

    ax.axis('off')

    return plt