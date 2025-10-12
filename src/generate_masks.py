#!/usr/bin/env python3
"""
Script:          generate_mask.py
Purpose:         Generates tissue masks, metadata, and QC plots
Author:          Sophia Li
Affiliation:     Campbell Lab
Date:            010-03-2025
"""

import argparse
from utils.preprocessing import preprocess_image
from utils.io_utils import (
    load_input_paths, 
    load_config, 
    load_panel, 
    load_image, 
    load_mask_generation_panel,
    save_tissue_mask,
    save_mask_metadata,
    save_mask_qc,
    generate_wsi_id_mapping
)
from utils.mask_utils import (
    determine_gmm_tissue_threshold, 
    determine_otsu_tissue_threshold, 
    generate_tissue_mask,
    generate_mask_qc_plot
)


def parse_arguments():
    # Parses and returns command-line arguments

    parser = argparse.ArgumentParser(
        description = "Generate tissue masks from WS-IMC images."
    )

    parser.add_argument("--config", type = str, required = True,
                        help = "Path to configuration file")
    parser.add_argument("--threshold-method", type = str, 
                        choices = ['otsu', 'gmm'], default = 'otsu',
                        help = "Thresholding method to determine a tissue threshold for the tissue mask: 'otsu', or 'gmm'")
    parser.add_argument("--remove-small-objects", type = bool, default = True,
                        help = "Remove small objects from the tissue mask")
    parser.add_argument("--fill-small-holes", type = bool, default = True,
                        help = "Fill small holes in the tissue mask")
    parser.add_argument("--save-metadata", type = bool, default = True,
                        help = "Save metadata associated with mask generation")
    parser.add_argument("--save-qc", type = bool, default = True, 
                        help = "Save quality control plots visualizing tissue mask coverage")
    parser.add_argument('--preprocess-images', type = bool, default = True,
                        help = "Preprocess the images in the input folder")
    return parser.parse_args()


def main():

    # Initialize the environment and parse user configurations
    args = parse_arguments()
    config = load_config(args.config)
    img_files, file_type = load_input_paths(config)
    id_mapper = generate_wsi_id_mapping(img_files, config)
    panel = load_panel(img_files[0], file_type, config)

    # Fetch mask generation markers if provided
    mask_panel = load_mask_generation_panel(config, panel)
    mask_metals = mask_panel['canonical_metal_tag'].tolist()

    # Generate a tissue mask for each input image
    for path in img_files:

        # Map the filepath to an internal WSI ID
        img_id = id_mapper[str(path)]

        # Load the image and filter for markers used for mask generation
        img = load_image(path, file_type, panel)
        composite = img.sel(channel = img['metal_tag'].isin(mask_metals))

        # Preprocess the image if toggled
        if args.preprocess_images:
            composite = preprocess_image(composite, config)

        # Dynamically generate a tissue threshold for the image: GMM or Otsu's
        if args.threshold_method == 'otsu':
            threshold, threshold_metadata  = determine_otsu_tissue_threshold(
                composite, config
            )
        else:
            threshold, threshold_metadata = determine_gmm_tissue_threshold(
                composite, config
            )

        # Generate a tissue mask using the threshold
        mask, mask_metadata = generate_tissue_mask(
            composite, threshold, args.remove_small_objects, 
            args.fill_small_holes, config
        )

        # Save the mask as the image's filename with '_mask.tiff' suffix
        save_tissue_mask(mask, img_id, config)

        # Save the mask generation metadata if toggled
        if args.save_metadata:
            metadata = {**threshold_metadata, **mask_metadata}
            save_mask_metadata(metadata, img_id, config)

        # Generate and save a quality control (QC) plot if toggled
        if args.save_qc:
            qc_plot = generate_mask_qc_plot(
                mask, composite, img_id, mask_panel,
                threshold_metadata, mask_metadata, config
            )
            save_mask_qc(qc_plot, img_id, config)


if __name__ == '__main__':
    main()