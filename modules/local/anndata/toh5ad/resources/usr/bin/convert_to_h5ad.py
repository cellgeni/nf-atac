#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import scanpy as sc


# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    stream=sys.stdout,  # Direct output to stdout instead of a file
)


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Script validates sample and annotation tables and splits annotation table into separate celltypes"
    )
    parser.add_argument(
        "--input",
        type=str,
        metavar="<path>",
        help="Specify a path to input file",
        required=True,
    )
    parser.add_argument(
        "--sample_id",
        type=str,
        metavar="<str>",
        default=None,
        help="Specify sample name for the file",
    )
    parser.add_argument(
        "--output",
        metavar="<path>",
        type=str,
        help="Specify a path to output .h5ad file",
        default=10,
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        metavar="<str>",
        default=None,
        help="Specify sample delimiter",
    )

    return parser


def check_10x_mtx_files(directory: str) -> bool:
    """
    Check if the directory contains the required 10x mtx files
    """
    required_files = ["matrix.mtx", "barcodes.tsv", "features.tsv"]
    for file in required_files:
        filepath = os.path.join(directory, file)
        if not os.path.isfile(filepath) and not os.path.isfile(filepath + ".gz"):
            logging.error(f"Missing required file: {file}")
            return False
    return True


def main() -> None:
    # Parse arguments
    parser = init_parser()
    args = parser.parse_args()

    # Check if input is directory
    if os.path.isdir(args.input):
        # return error if any of the required files are missing
        if not check_10x_mtx_files(args.input):
            raise FileNotFoundError(
                "The specified directory does not contain the required files: matrix.mtx, barcodes.tsv, and features.tsv"
            )

        # load 10x mtx file
        logging.info("Loading 10x mtx file to AnnData object")
        adata = sc.read_10x_mtx(
            args.input,
            var_names="gene_symbols",
            gex_only=True,
        )
    else:
        # get file extension
        _, extension = os.path.splitext(args.input)

        # read file based on extension using case match
        match extension:
            case ".h5":
                logging.info("Loading .h5 file to AnnData object")
                adata = sc.read_10x_h5(args.input, gex_only=True)
            case ".mtx":
                logging.info("Loading .mtx file to AnnData object")
                adata = sc.read_mtx(args.input).T
            case ".zarr":
                logging.info("Loading .zarr file to AnnData object")
                raise NotImplementedError(
                    "Loading .zarr file is not implemented yet. Please provide a path to .h5 file or .mtx file."
                )
            case _:
                raise ValueError(
                    "Unsupported file format. Please provide a path to .h5 file, .mtx file or .mtx file directory."
                )

    # Make var_names unique
    adata.var_names_make_unique()

    # Add sample name to obs
    logging.info("Adding sample name to obs")
    adata.obs["sample"] = args.sample_id

    # Add delimiter to obs index if specified
    if args.delimiter:
        logging.info("Adding delimiter to obs index")
        adata.obs["barcode"] = adata.obs.index
        adata.obs.index = adata.obs["barcode"] + args.delimiter + adata.obs["sample"]
        adata.obs.index.name = "barcode_sample"

    # Save adata abject
    logging.info("Saving AnnData object to .h5ad file")
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
