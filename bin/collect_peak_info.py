#!/usr/bin/env python3

import os
import argparse
from pandas import DataFrame


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Derives consensus peaks from narrow peaks using the TGCA iterative peak filtering approach"
    )
    parser.add_argument(
        "--publish_dir",
        type=str,
        metavar="<dir>",
        help="Path to publishing directory",
    )
    parser.add_argument(
        "--celltype_names",
        metavar="<val>",
        type=str,
        nargs="+",
        help="Specify all celltype names to be listed in the metadata file",
    )
    parser.add_argument(
        "--fragment_counts",
        metavar="<val>",
        type=int,
        nargs="+",
        help="Specify fragment counts to the specified celltypes",
    )
    parser.add_argument(
        "--large_peak_counts",
        metavar="<val>",
        type=int,
        nargs="+",
        help="Specify large peak counts to the specified celltypes",
    )
    parser.add_argument(
        "--all_peak_counts",
        metavar="<val>",
        type=int,
        nargs="+",
        help="Specify peak counts to the specified celltypes",
    )
    parser.add_argument(
        "--narrowPeaks",
        metavar="<file>",
        type=str,
        nargs="+",
        help="Specify names for .narrowPeak files to the specified celltypes",
    )
    parser.add_argument(
        "--output",
        metavar="<file>",
        type=str,
        help="Specify a name to the output file",
        default="pseudobulk_peaks.tsv",
    )
    return parser


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # read inputs to Dict file to DataFrame
    metadata_dict = {
        "celltype_names": args.celltype_names,
        "fragment_counts": args.fragment_counts,
        "large_peak_counts": args.large_peak_counts,
        "all_peak_counts": args.all_peak_counts,
        "filepath": [os.path.join(args.publish_dir, path) for path in args.narrowPeaks],
    }

    # convert to DataFrame
    metadata_df = DataFrame(metadata_dict)

    # save to file
    metadata_df.to_csv(args.output, index=False, sep='\t')


if __name__ == "__main__":
    main()
