#!/usr/bin/env python
import argparse
import os
import glob
from typing import Dict
from make_pseudobulk import read_chromsizes
import logging
from colored_logger import setup_logging
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pyranges import PyRanges
from pandas import read_csv
from pandas.errors import DtypeWarning
import warnings

# ignore pandas DtypeWarning warnings
warnings.simplefilter(action="ignore", category=DtypeWarning)


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
        "--sample_id", type=str, metavar="<val>", help="Sample identificator"
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
    return parser


def path_to_celltype(filepath: str) -> str:
    """
    Get celltype name from .narrowPeak file path
    filepath (str): a path to the .narrowPeak file
    """
    filename = os.path.basename(filepath)
    celltype = filename.replace("_peaks.narrowPeak", "")
    return celltype


def check_if_empty(filepath: str, sample_id: str, skip_empty_peaks: bool) -> bool:
    """
    Checks if .narrowPeak file is empty
    filepath (str): A path to the .narrowPeak file
    sample_id (str): Sample identifier
    skip_empty_peaks (bool): If True skips celltypes with no peaks found
    """
    # check if file is empty
    file_size = os.stat(filepath).st_size

    if file_size == 0 and skip_empty_peaks:
        logging.warning(
            f"{sample_id} has no peaks for {path_to_celltype(filepath)}. Skipping"
        )
        return True
    elif file_size == 0 and not skip_empty_peaks:
        raise ValueError(
            f"{sample_id} has no peaks for {path_to_celltype(filepath)}, exiting. Set skip_empty_peaks to True to skip empty peaks."
        )
    return False


def read_narrow_peak(filepath: str) -> PyRanges:
    """
    Reads .narrowPeak file to PyRanges objects
    filepath (str): A path to the .narrowPeak file
    sample_id (str): Sample identifier
    skip_empty_peaks (bool): If True skips celltypes with no peaks found
    """
    # read narrow peaks to DataFrame
    narrow_peak = read_csv(filepath, sep="\t", header=None)
    narrow_peak.columns = NARROW_PEAK_COLUMNS

    # convert to PyRanges object
    narrow_peak_pr = PyRanges(narrow_peak)
    return narrow_peak_pr


def main():
    # set up logger
    setup_logging()

    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # read chromsizes file to DataFrame
    chromsizes = read_chromsizes(args.chromsizes)

    # read .narrowpeaks files to pr.PyRanges objects
    narrow_peak_files = glob.glob(os.path.join(args.narrow_peaks, "*.narrowPeak"))
    narrow_peak_dict = {
        path_to_celltype(filepath): read_narrow_peak(filepath)
        for filepath in narrow_peak_files
        if not check_if_empty(filepath, args.sample_id, args.skip_empty_peaks)
    }

    # infer consensus peaks
    consensus_peaks = get_consensus_peaks(
        narrow_peaks_dict=narrow_peak_dict,
        peak_half_width=args.peak_half_width,
        chromsizes=chromsizes,
        path_to_blacklist=args.blacklist,
    )

    # save consensus peaks to .bed
    consensus_peaks.to_bed(
        path=f"{args.sample_id}_consensus.bed",
        keep=True,
        compression="infer",
        chain=False,
    )


if __name__ == "__main__":
    main()
