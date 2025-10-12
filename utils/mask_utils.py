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
import matplotlib.gridspec as gridspec

def determine_otsu_tissue_threshold(composite, config):
    # Determines and returns a global intensity threshold to separate tissue 
    # from background using Otsu's method

    # Extract the values of the xarray composite
    composite = composite.values

    # Fetch configurations
    min_tissue_threshold = (config
                            .get('tissue_mask', {})
                            .get('min_tissue_threshold', 0))
    
    # Generate the threshold using Otsu's method and enforce a minimum
    threshold = threshold_otsu(composite.flatten())
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
                            .get('tissue_mask', {})
                            .get('min_tissue_threshold', 0))
    seed = config.get('seed', 42)

    # Aggregate across channels using the median and flatten to pixels
    composite = np.max(composite, axis = 0)
    pixels = composite.values.flatten().reshape(-1, 1)

    # # Cip extreme outliers that may skew background Gaussian
    # pixels = np.clip(
    #     pixels, np.percentile(pixels, 1), np.percentile(pixels, 99)
    # )
    
    # Ignore background to avoid skewing towards lower values
    pixels = pixels[pixels > 0.01].reshape(-1, 1)

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
        #threshold = max(threshold, np.percentile(pixels, 50))
        metadata = {
            'method': 'gmm',
            'threshold': threshold,
            'min_threshold': min_tissue_threshold,
            'gmm_means_1': means[0],
            'gmm_means_2': means[1],
            'gmm_stds_1': stds[0],
            'gmm_stds_2': stds[1],
            'gmm_weights_1': weights[0],
            'gmm_weights_2': weights[1]
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
    mask = (composite > threshold).any(axis = 0).values
    orig_area_px = np.sum(mask)

    original_mask = mask.copy()

    # Remove small objects if toggled, measuring them
    if remove_objects:
        mask = remove_small_objects(mask, object_threshold)
        removed_objects = original_mask & ~mask
        total_removed_area = np.sum(removed_objects)
    else:
        total_removed_area = 0

    # Fill small holes if toggled, measuring them
    if fill_holes:
        mask = remove_small_holes(mask, hole_threshold)
        filled_holes = mask & ~original_mask
        total_filled_area = np.sum(filled_holes)
    else:
        total_filled_area = 0

    # Record mask metadata
    final_area_px = np.sum(mask)
    img_area_px = composite.shape[1] * composite.shape[2]

    metadata = {
        'raw_mask_area_px': int(orig_area_px),
        'clean_mask_area_px': int(final_area_px),
        'full_img_area_px': int(img_area_px),
        "removed_objects_area_px": int(total_removed_area),
        "filled_holes_area_px": int(total_filled_area),
        "mask_coverage_percent": 100 * final_area_px / img_area_px,
        "cleaning_area_change_px": int(final_area_px - orig_area_px),
        "cleaning_area_change_percent": ((final_area_px - orig_area_px) 
                                / orig_area_px * 100 
                                if orig_area_px > 0 else 0)
    }

    return mask, metadata


def generate_mask_qc_plot(mask, composite, image_name, mask_panel,
                          threshold_metadata, mask_metadata, config):
    # Generates and returns a quality control (QC) plot overlaying the tissue
    # mask upon the composite image

    # Filter the composite for the RGB markers specified in the configurations
    rgb_markers = config.get('tissue_mask', {}).get('rgb_markers')
    rgb_panel = mask_panel[mask_panel['canonical_marker'].isin(rgb_markers)]
    rgb_metal_tags = rgb_panel['canonical_metal_tag'].tolist()

    composite = composite.sel(
        channel = np.isin(composite['metal_tag'], rgb_metal_tags)
    ).values
    composite = np.transpose(composite, (1, 2, 0))

    # Plot and return the figure: image on the left, overlaid on the right
    fig = plt.figure(figsize = (14, 7))
    gs = gridspec.GridSpec(1, 2, figure = fig, width_ratios = [1, 1], 
                           wspace = 0.01, left = 0.1, right = 0.9, top = 0.9)  
    
    ax0, ax1 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # Left: composite image only
    ax0.imshow(composite)
    ax0.axis('off')
    ax0.set_position([0.152, 0.1, 0.38, 0.78])

    # Right: composite with tissue mask overlaid in transparent red
    overlay = np.zeros((*mask.shape, 4))
    overlay[mask > 0] = [1, 0, 0, 0.4]
    ax1.imshow(composite)
    ax1.imshow(overlay)
    ax1.axis('off')
    ax1.set_position([0.4855, 0.1, 0.38, 0.78])

    # Set titles and subplot titles
    fig.text(0.5, 0.925, f"Tissue Mask QC - {image_name}", 
             ha = 'center', va = 'top', fontsize = 16,
             transform = fig.transFigure)
    fig.text(ax0.get_position().x0 + ax0.get_position().width / 2, 0.065, 
             "Composite", ha = 'center', fontsize = 12)
    fig.text(ax1.get_position().x0 + ax1.get_position().width / 2, 0.065, 
             "Composite with Mask Overlay", ha = 'center', fontsize = 12)
    
    # Build metadata and RGB marker captions
    threshold_caption = (
        f"Method: {threshold_metadata['method']}    "
        f"Threshold: {threshold_metadata['threshold']:.4f}  "
        f"Coverage: {mask_metadata['mask_coverage_percent']:.2f}%"
    )
    rgb_caption = (
        f"Red: {rgb_markers[0]} ({rgb_metal_tags[0]})   "
        f"Green: {rgb_markers[1]} ({rgb_metal_tags[1]})     "
        f"Blue: {rgb_markers[2]} ({rgb_metal_tags[2]})"
    )

    fig.text(0.32, 0.03, threshold_caption, ha = 'center',
             fontsize = 9, color = 'black')
    fig.text(0.65, 0.03, rgb_caption, ha = 'center',
             fontsize = 9, color = 'black')

    return fig