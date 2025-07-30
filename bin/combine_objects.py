#!/usr/bin/env python3

import os
import pickle
import glob
import argparse
import anndata as ad
from pycisTopic.cistopic_class import merge, CistopicObject


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Creates cisTopic object from fragments, consensus peaks and quality control files"
    )
    parser.add_argument(
        "--cistopic",
        metavar="<dir>",
        type=str,
        nargs="+",
        help="Specify a path to the directory with cistopic objects",
    )
    parser.add_argument(
        "--anndata",
        metavar="<dir>",
        type=str,
        nargs="+",
        help="Specify a path to the directory with anndata objects",
    )
    parser.add_argument(
        "--combined_cistopic",
        metavar="<file>",
        type=str,
        help="Specify a name to the file with combined cistopic object",
    )
    parser.add_argument(
        "--combined_anndata",
        metavar="<file>",
        type=str,
        help="Specify a name to the file with combined cistopic object",
    )
    parser.add_argument(
        "--output_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to output directory",
    )
    return parser


def read_cistopic_object(file: str) -> CistopicObject:
    """
    Read cistopic object from the file
    """
    with open(file, "rb") as f:
        cistopic_obj = pickle.load(f)
    return cistopic_obj


def main():
    """
    Main function of the script
    """
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # combine cistopic objects
    cistopic_objects = [read_cistopic_object(file) for file in args.cistopic]
    cistopic_combined = merge(cistopic_objects)

    # write cistopic object to the file
    with open(os.path.join(args.output_dir, args.combined_cistopic), "wb") as f:
        pickle.dump(cistopic_combined, f)

    # delete cistopic objects from memory
    del cistopic_objects
    del cistopic_combined

    # combine anndata objects and write to disk
    anndata_objects = [ad.read_h5ad(file) for file in args.anndata]
    adata_combined = ad.concat(
        adatas=anndata_objects,
        join="outer",
    )

    adata_combined.write_h5ad(
        os.path.join(args.output_dir, args.combined_anndata),
    )


if __name__ == "__main__":
    main()
