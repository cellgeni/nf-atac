#!/usr/bin/env python3

import argparse
import anndata as ad


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Creates cisTopic object from fragments, consensus peaks and quality control files"
    )
    parser.add_argument(
        "--anndata",
        metavar="<file>",
        type=str,
        nargs="+",
        help="Specify a path to the .h5ad files",
    )
    parser.add_argument(
        "--output",
        metavar="<file>",
        type=str,
        help="Specify a name to the file with combined AnnData object",
    )
    return parser


def main():
    """
    Main function of the script
    """
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # combine anndata objects
    anndata_objects = [ad.read_h5ad(file) for file in args.anndata]
    anndata = ad.concat(anndata_objects, join="outer", merge="same")

    # write anndata object to the file
    anndata.write_h5ad(args.output)


if __name__ == "__main__":
    main()