#!/usr/bin/env python3

import os
import glob
import logging
import warnings
import argparse
from typing import Dict
from pandas import read_csv
from pyranges import PyRanges
from colored_logger import setup_logging
from make_pseudobulk import read_chromsizes
from pycisTopic.iterative_peak_calling import get_consensus_peaks


NARROW_PEAK_COLUMNS = [
    "Chromosome",
    "Start",
    "End",
    "Name",
    "Score",
    "Strand",
    "FC_summit",
    "-log10_pval",
    "-log10_qval",
    "Summit",
]


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Derives consensus peaks from narrow peaks using the TGCA iterative peak filtering approach"
    )
    parser.add_argument(
        "--narrow_peaks",
        type=str,
        metavar="<dir>",
        help="Path to dirrectory with narrow peaks files in the format celltype_peaks.narrowPeak",
    )
    parser.add_argument(
        "--chromsizes",
        metavar="<file>",
        type=str,
        help="Specify a path to the file with chromosome lengths from the UCSC databases",
    )
    parser.add_argument(
        "--blacklist",
        metavar="<file>",
        type=str,
        help="Specify a path to bed file containing blacklist regions (Amemiya et al., 2019)",
    )
    parser.add_argument(
        "--consensus",
        metavar="<file>",
        type=str,
        help="Specify an output file name",
    )
    parser.add_argument(
        "--peak_half_width",
        metavar="<val>",
        type=int,
        help="Number of base pairs that each summit will be extended in each direction",
        default=250,
    )
    parser.add_argument(
        "--skip_empty_peaks",
        help="If specified skips celltypes with no peaks found",
        action="store_true",
    )
    parser.add_argument(
        "--logfile",
        metavar="<file>",
        type=str,
        help="Specify a log file name",
        default="inferconsensus.log",
    )
    return parser


def path_to_celltype(filepath: str) -> str:
    """
    Get celltype name from .narrowPeak file path
    filepath (str): a path to the .narrowPeak file
    """
    filename = os.path.basename(filepath)
    celltype = filename.replace("_peaks.narrowPeak", "")
    return celltype


def check_if_empty(filepath: str, skip_empty_peaks: bool) -> bool:
    """
    Checks if .narrowPeak file is empty
    filepath (str): A path to the .narrowPeak file
    skip_empty_peaks (bool): If True skips celltypes with no peaks found
    """
    # check if file is empty
    file_size = os.stat(filepath).st_size

    if file_size == 0 and skip_empty_peaks:
        logging.warning(f"No peaks found for {path_to_celltype(filepath)}. Skipping")
        return True
    elif file_size == 0 and not skip_empty_peaks:
        raise ValueError(
            f"No peaks found for {path_to_celltype(filepath)}, exiting. Set skip_empty_peaks to True to skip empty peaks."
        )
    return False


def read_narrow_peak(filepath: str) -> PyRanges:
    """
    Reads .narrowPeak file to PyRanges objects
    filepath (str): A path to the .narrowPeak file
    """
    # read narrow peaks to DataFrame
    narrow_peak = read_csv(filepath, sep="\t", header=None)
    narrow_peak.columns = NARROW_PEAK_COLUMNS

    # convert to PyRanges object
    narrow_peak_pr = PyRanges(narrow_peak)
    return narrow_peak_pr


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # set up logger
    setup_logging(args.logfile)

    # read chromsizes file to DataFrame
    chromsizes = read_chromsizes(args.chromsizes)

    # read .narrowpeaks files to pr.PyRanges objects
    narrow_peak_files = glob.glob(os.path.join(args.narrow_peaks, "*.narrowPeak"))
    narrow_peak_dict = {
        path_to_celltype(filepath): read_narrow_peak(filepath)
        for filepath in narrow_peak_files
        if not check_if_empty(filepath, args.skip_empty_peaks)
    }

    # infer consensus peaks
    with warnings.catch_warnings(action="ignore"):
        consensus_peaks = get_consensus_peaks(
            narrow_peaks_dict=narrow_peak_dict,
            peak_half_width=args.peak_half_width,
            chromsizes=chromsizes,
            path_to_blacklist=args.blacklist,
        )

    # save consensus peaks to .bed
    consensus_peaks.to_bed(
        path=args.consensus,
        keep=True,
        compression="infer",
        chain=False,
    )


if __name__ == "__main__":
    main()
