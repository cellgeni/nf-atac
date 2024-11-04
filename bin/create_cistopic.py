#!/usr/bin/env python
import os
import pickle
import argparse
from typing import List
from polars import read_parquet
from numpy import ndarray
from pandas import DataFrame
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Creates cisTopic object from fragments, consensus peaks and quality control files"
    )
    parser.add_argument(
        "--sample_id", type=str, metavar="<val>", help="Sample identificator"
    )
    parser.add_argument(
        "--fragments",
        metavar="<file>",
        type=str,
        help="Specify a path to the fragments.tsv.gz file (fragments.tsv.gz should be in the same directory)",
    )
    parser.add_argument(
        "--consensus",
        metavar="<file>",
        type=str,
        help="Specify a path to the file with consensus peaks",
    )
    parser.add_argument(
        "--blacklist",
        metavar="<file>",
        type=str,
        help="Specify a path to bed file containing blacklist regions (Amemiya et al., 2019)",
    )
    parser.add_argument(
        "--qc_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to the directory with qualtiry control results",
    )
    parser.add_argument(
        "--cpus", metavar="<num>", type=int, help="Specify a number of cpu cores to use"
    )
    parser.add_argument(
        "--unique_fragments_threshold",
        metavar="<val>",
        type=int,
        help="Threshold for number of unique fragments in peaks",
    )
    parser.add_argument(
        "--tss_enrichment_threshold",
        metavar="<val>",
        type=float,
        help="Threshold for TSS enrichment score",
    )
    parser.add_argument(
        "--frip_threshold",
        metavar="<val>",
        type=float,
        help="Threshold for fraction of reads in peaks (FRiP). If not defined the threshold will be set to 0",
        default=0,
    )
    parser.add_argument(
        "--use_automatic_thresholds",
        help="Use automatic thresholds for unique fragments in peaks and TSS enrichment score as calculated by Otsu's method",
        action="store_true",
    )
    parser.add_argument(
        "--split_pattern",
        metavar="<val>",
        type=str,
        help="Pattern to split cell barcode from sample id. Default: '___'",
        default="___",
    )
    return parser


def read_metrics(qc_dir: str, sample_id: str, barcodes: ndarray) -> DataFrame:
    """
    Read QC metrics from <sample_id>.fragments_stats_per_cb.parquet file
    qc_dir (str): directory with QC results
    sample_id (str): sample identifier
    barcodes (List[str]): a list of barcodes that passed filtering
    """
    # get file's path
    fragments_stats_file = os.path.join(
        qc_dir, f"{sample_id}.fragments_stats_per_cb.parquet"
    )
    if not os.path.exists(fragments_stats_file):
        raise ValueError('No file with path "{fragments_stats_file}" was found')
    # if there is only 1 passing barcode it is returned as 0d array
    # as result of numpy squeeze function see
    # https://github.com/aertslab/pycisTopic/blob/787ce422a37f5975b0ebb9e7b19eeaed44847501/src/pycisTopic/qc.py#L140C40-L140C50
    # barcodes = barcodes if barcodes.shape else [barcodes[()]]
    # read data and filter barcodes
    fragments_stats = read_parquet(fragments_stats_file).to_pandas()
    fragments_stats = fragments_stats.set_index("CB")
    fragments_stats = fragments_stats.loc[barcodes].copy()
    return fragments_stats


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # filter low quality cells
    barcodes, thresholds = get_barcodes_passing_qc_for_sample(
        sample_id=args.sample_id,
        pycistopic_qc_output_dir=args.qc_dir,
        unique_fragments_threshold=args.unique_fragments_threshold,
        tss_enrichment_threshold=args.tss_enrichment_threshold,
        frip_threshold=args.frip_threshold,
        use_automatic_thresholds=args.use_automatic_thresholds,
    )

    # read fragments stats
    fragments_stats = read_metrics(args.qc_dir, args.sample_id, barcodes)

    # create a cistopic object
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments=args.fragments,
        path_to_regions=args.consensus,
        path_to_blacklist=args.blacklist,
        metrics=fragments_stats,
        valid_bc=barcodes,
        n_cpu=args.cpus,
        project=args.sample_id,
        split_pattern=args.split_pattern,
    )

    # save to pickle file
    with open(f"{args.sample_id}_cistopic_obj.pkl", "wb") as file:
        pickle.dump(cistopic_obj, file)


if __name__ == "__main__":
    main()
