#!/usr/bin/env python3

import argparse

from utils.io_utils import (
    load_input_paths, 
    load_config, load_panel, 
    load_image, 
    load_mask_generation_panel,
    save_tissue_mask,
    save_mask_metadata,
    save_mask_qc
)
from utils.mask_utils import (
    determine_gmm_tissue_threshold, 
    determine_otsu_tissue_threshold, 
    generate_tissue_mask,
    generate_mask_qc_plot
)
from utils.preprocessing import preprocess_image


from skimage import morphology, measure, exposure
from skimage.measure import find_contours


def parse_arguments():
    # Parses and returns command-line arguments

    parser = argparse.ArgumentParser(
        description = "Generate tissue masks from WS-IMC images."
    )
    # parser.add_argument("--input", type = str, required = True, 
    #                     help = "Path to input folder")
    # parser.add_argument("--output", type = str, required = True,
    #                     help = "Path to output folder")
    parser.add_argument("--config", type = str, default = 'config.yaml',
                        help = "Path to configuration file")
    # parser.add_argument("--panel", type = str, required = True,
    #                     help = "Path to panel")
    parser.add_argument("--threshold-method", type = str, 
                        choices = ['otsu', 'gmm'], default = 'otsu',
                        help = "Thresholding method to determine a tissue threshold for the tissue mask: 'otsu', or 'gmm'")
    parser.add_argument("--remove-small-objects", type = bool, default = True,
                        help = "Remove small objects from the tissue mask")
    parser.add_argument("--fill-small-holes", type = bool, default = True,
                        help = "Fill small holes in the tissue mask")
    parser.add_argument("--save_metadata", type = bool, default = True,
                        help = "Save metadata associated with mask generation")
    parser.add_argument("--save_qc", type = bool, default = True, 
                        help = "Save quality control plots visualizing tissue mask coverage")
    parser.add_argument('--preprocess_images', type = bool, default = True,
                        help = "Preprocess the images in the input folder")
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
        image_name = path.name

        # Load the image and preprocess if toggled
        image = load_image(path)
        if args.preprocess_image:
            composite = preprocess_image(image, config)

        # Filter for markers used for mask generation after preprocessing
        composite = image.sel(metal_tag = mask_metals)

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
        save_tissue_mask(mask, image_name, config)

        # Save the mask generation metadata if toggled
        if args.save_metadata:
            metadata = {**threshold_metadata, **mask_metadata}
            save_mask_metadata(metadata, image_name, config)

        # Generate and save a quality control (QC) plot if toggled
        if args.save_qc:
            qc_plot = generate_mask_qc_plot(
                mask, composite, image_name, threshold_metadata, 
                mask_metadata, config
            )
            save_mask_qc(qc_plot, image_name, config)


