#!/usr/bin/env python3

import argparse
import os
import sys
import logging
from typing import List, Dict
from colored_logger import setup_logging
from tempfile import TemporaryDirectory
from pandas import read_csv, read_table, DataFrame, Series
from pandas.errors import InvalidColumnName, DataError
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Create pseudobulks as bed and bigwig from single cell fragments file given a barcode annotation"
    )
    parser.add_argument(
        "--sample_id",
        type=str,
        nargs="+",
        help="Specify sample ids in celltype annotation file you want to process",
    )
    parser.add_argument(
        "--fragments",
        metavar="<file>",
        type=str,
        nargs="+",
        help="Specify a paths to the fragments.tsv.gz files in the same order as sample ids (fragments.tsv.gz should be in the same directory)",
    )
    parser.add_argument(
        "--celltype_annotation",
        metavar="<file>",
        type=str,
        help="Specify a path to the cell-type annotation file",
    )
    parser.add_argument(
        "--fragments_celltype_x_sample",
        metavar="<file>",
        type=str,
        help="Specify a path to the fragments_celltype_x_sample.csv file",
    )
    parser.add_argument(
        "--chromsizes",
        metavar="<file>",
        type=str,
        help="Specify a path to the file with chromosome lengths from the UCSC databases",
    )
    parser.add_argument(
        "--output_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to the output directory",
    )
    parser.add_argument(
        "--bigwig_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to the bigwig output directory",
        default="bigwig",
    )
    parser.add_argument(
        "--celltype_col",
        metavar="<val>",
        type=str,
        help="Specify a name for celltype column in annotation file",
        default="celltype",
    )
    parser.add_argument(
        "--sample_id_col",
        metavar="<val>",
        type=str,
        help="Specify a name for sample_id column in annotation file",
        default="sample_id",
    )
    parser.add_argument(
        "--barcode_col",
        metavar="<val>",
        type=str,
        help="Specify a name for barcode column in annotation file and per_barcode_metrics file",
        default="barcode",
    )
    parser.add_argument(
        "--cpus", metavar="<num>", type=int, help="Specify a number of cpu cores to use"
    )
    parser.add_argument(
        "--logfile",
        metavar="<file>",
        type=str,
        help="Specify a log file name",
        default="pseudobulk.log",
    )
    return parser


def read_chromsizes(chromsizes_file: str) -> DataFrame:
    """
    Read chromsizes file to pandas DataFrame and add "Start" column
    Args:
        chromsizes_file (str): chromsizes_file (str): path to the file with chromosome lengths from the UCSC databases

    Returns:
        DataFrame: contains chromsizes
    """
    chromsizes = read_table(chromsizes_file, header=None, names=["Chromosome", "End"])
    chromsizes.insert(1, "Start", 0)
    return chromsizes


def check_missing_samples(
    celltypes: DataFrame, sample_list: List[str], sample_id_col: str, celltype_name: str
):
    """
    Check if there are samples with no barcodes available
    Args:
        celltypes (DataFrame): a dataframe with barcode/celltype annotation
        sample_list (List[str]): a list of samples specified in sample table
        sample_id_col (str): a sample_id column in annotation file
        celltype_name (str): a celltype name
    """
    samples_in_annotation = celltypes[sample_id_col].unique().tolist()
    difference = set(sample_list).difference(samples_in_annotation)
    if difference:
        sample_string = ",".join(difference)
        logging.warning(
            f"Following samples have zero barcodes for '{celltype_name}' celltype: {sample_string}"
        )


def get_fragments_per_sample(
    filepath: str,
    samples_in_annotation: List[str],
    celltype_name: str,
    celltype_col: str,
) -> Series:
    """
    Reads fragments_celltype_x_sample.csv file
    Args:
        filepath (str): a path to fragments_celltype_x_sample.csv file
        samples_in_annotation (List[str]): a list of sample names that we want to read from fragments_per_sample.csv file
        celltype_name (str): a name of the celltype that we want to subsample
        celltype_col (str): a name of celltype column in annotation file

    Returns:
        Series: sample fragments for `celltype_name`
    """
    columns_to_use = [celltype_col] + samples_in_annotation
    fragments_celltype_x_sample = read_csv(
        filepath, index_col=0, usecols=columns_to_use
    )
    fragments_per_sample = fragments_celltype_x_sample.loc[celltype_name].copy()
    return fragments_per_sample


def get_zero_fragment_samples(
    fragments_celltype_x_sample_path: str,
    celltypes: DataFrame,
    celltype_name: str,
    sample_id_col: str,
    celltype_col: str,
) -> List[str]:
    """
    Get a sample names with zero fragments found for the `celltype_name`
    Args:
        fragments_celltype_x_sample_path (str): a path to fragments_celltype_x_sample_path.csv file
        celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
        celltype_name (str):  a name of the celltype in annotation file
        sample_id_col (str): a name of the column with sample_id labels in celltypes DataFrame
        celltype_col (str): a name of celltype column in annotation file

    Returns:
        List[str]: a sample list with zero fragments found for the `celltype_name`
    """
    # get samples from annotation
    samples_in_annotation = celltypes[sample_id_col].unique().tolist()
    # read fragments_celltype_x_sample file
    fragments_per_sample = get_fragments_per_sample(
        fragments_celltype_x_sample_path,
        samples_in_annotation,
        celltype_name,
        celltype_col,
    )
    # get samples with zero fragments
    zero_fragment_samples = fragments_per_sample[
        fragments_per_sample == 0
    ].index.tolist()
    return zero_fragment_samples


def drop_empty_samples(
    sample_list: List[str], celltypes: DataFrame, celltype_name: str, sample_id_col: str
) -> DataFrame:
    """
    Drops sample with zero fragments from celltype DataFrame
    Args:
        sample_list (List[str]): a list of samples that we want to drop
        celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
        celltype_name (str): a name of the celltype in annotation file
        sample_id_col (str): a name of the column with sample_id labels in celltypes DataFrame

    Returns:
        DataFrame: a celltype DataFrame without samples from sample_list
    """
    if sample_list:
        sample_string = ",".join(sample_list)
        # log this step
        logging.info(
            f"Dropping {sample_string} from the from {celltype_name} annotation file as no fragments is found"
        )
        # drop the celltypes
        zero_fragment_cells = celltypes[
            celltypes[sample_id_col].isin(sample_list)
        ].index
        celltypes_filtered = celltypes.drop(zero_fragment_cells)
        return celltypes_filtered
    return celltypes


def filter_empty_fragments(
    celltypes: DataFrame,
    fragments_celltype_x_sample: str,
    celltype_name: str,
    sample_id_col: str,
    celltype_col: str,
) -> DataFrame:
    """
    Filter samples with zero fragments found
    Args:
        celltypes (DataFrame): a dataframe with barcode/celltype annotation
        fragments_celltype_x_sample (str): a path to fragments_celltype_x_sample.csv file
        celltype_name (str): a name of the celltype in annotation file
        sample_id_col (str): a name of the column with sample_id labels in celltypes DataFrame
        celltype_col (str): a name of celltype column in annotation file

    Returns:
        DataFrame: filtered dataframe withou zero fragment samples
    """
    # get zero fragment samples
    zero_fragment_samples = get_zero_fragment_samples(
        fragments_celltype_x_sample,
        celltypes,
        celltype_name,
        sample_id_col,
        celltype_col,
    )
    # drop samples with zero fragments found
    filtered_celltypes = drop_empty_samples(
        zero_fragment_samples, celltypes, celltype_name, sample_id_col
    )
    return filtered_celltypes


def read_celltype_annotation(
    filepath: str,
    sample_list: List[str],
    fragments_celltype_x_sample: str,
    sample_id_col: str,
    celltype_col: str,
) -> DataFrame:
    """
    Read cell-type annotation file to pandas DataFrame
    Args:
        filepath (str): a path to celltype annotation
        sample_list (List[str]): a list of samples mentioned in sample table
        fragments_celltype_x_sample (str): a path to fragments_celltype_x_sample.csv file
        sample_id_col (str): a name of the column with sample_id labels in celltypes DataFrame
        celltype_col (str): a name of celltype column in annotation file

    Raises:
        SystemExit: if no fragments found for the celltype

    Returns:
        DataFrame: celltype annotation
    """
    # read chromsizes file
    celltype_name = filepath.replace(".csv", "").replace("_", " ")
    celltypes = read_csv(filepath)
    # check if all samples are present in annotation
    check_missing_samples(celltypes, sample_list, sample_id_col, celltype_name)
    # check if there are samples with zero fragments for this celltype
    filtered_celltypes = filter_empty_fragments(
        celltypes,
        fragments_celltype_x_sample,
        celltype_name,
        sample_id_col,
        celltype_col,
    )
    # check if DataFrame is empty
    if filtered_celltypes.empty:
        logging.info(f"There is no fragments found for {celltype_name} celltype")
        raise SystemExit
    # get sample list
    filtered_sample_list = filtered_celltypes[sample_id_col].unique().tolist()
    return filtered_celltypes, filtered_sample_list


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # set up logger
    setup_logging(args.logfile)

    # make a dictionaries for fragments
    fragments_dict = dict(zip(args.sample_id, args.fragments))

    # read files to DataFrame
    chromsizes = read_chromsizes(args.chromsizes)
    celltypes, sample_list = read_celltype_annotation(
        args.celltype_annotation,
        args.sample_id,
        args.fragments_celltype_x_sample,
        args.sample_id_col,
        args.celltype_col,
    )

    # subsample samples present in celltype DataFrame
    filtered_fragments_dict = {
        sample_id: fragments_dict[sample_id] for sample_id in sample_list
    }

    # make output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.bigwig_dir, exist_ok=True)

    # run pseudobulk generation
    with TemporaryDirectory(prefix="pb", dir="/tmp") as temp_dir:
        bw_paths, bed_paths = export_pseudobulk(
            input_data=celltypes,
            variable=args.celltype_col,
            sample_id_col=args.sample_id_col,
            chromsizes=chromsizes,
            bed_path=args.output_dir,
            bigwig_path=args.bigwig_dir,
            path_to_fragments=filtered_fragments_dict,
            n_cpu=args.cpus,
            normalize_bigwig=True,
            temp_dir=temp_dir,
        )


if __name__ == "__main__":
    main()
