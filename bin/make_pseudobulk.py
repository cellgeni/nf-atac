#!/usr/bin/env python
import argparse
import os
import sys
import logging
from typing import List
from colored_logger import setup_logging
from pandas import read_csv, read_table, DataFrame
from pandas.errors import InvalidColumnName, DataError
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

CELLTYPE_COLUMNS = {"sample_id", "barcode", "celltype"}


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
        help="Sample name in celltype annotation file",
    )
    parser.add_argument(
        "--fragments",
        metavar="<file>",
        type=str,
        help="Specify a path to the fragments.tsv.gz file (fragments.tsv.gz should be in the same directory)",
        default="atac_fragments.tsv.gz",
    )
    parser.add_argument(
        "--celltype_annotation",
        metavar="<file>",
        type=str,
        help="Specify a path to the cell-type annotation file",
    )
    parser.add_argument(
        "--chromsizes",
        metavar="<file>",
        type=str,
        help="Specify a path to the file with chromosome lengths from the UCSC databases",
    )
    parser.add_argument(
        "--barcode_metrics",
        metavar="<file>",
        type=str,
        help="Specify Cellrangers per_barcode_metrics.csv filename",
        default="per_barcode_metrics.csv",
    )
    parser.add_argument(
        "--output_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to the output directory",
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
        "--fragments_col",
        metavar="<val>",
        type=str,
        help="Specify a name for fragments column in 'per_barcode_metrics.csv' file",
        default="atac_fragments",
    )
    parser.add_argument(
        "--dropna",
        help="If specified drops NaN values from celltype annotation",
        action="store_true",
    )
    parser.add_argument(
        "--skip_empty_fragments",
        help="If specified skips celltypes with no fragments found",
        action="store_true",
    )
    parser.add_argument(
        "--cpus", metavar="<num>", type=int, help="Specify a number of cpu cores to use"
    )
    return parser


def read_chromsizes(chromsizes_file: str) -> DataFrame:
    """
    Read chromsizes file to pandas DataFrame and add "Start" column
    chromsizes_file (str): path to the file with chromosome lengths from the UCSC databases
    """
    chromsizes = read_table(chromsizes_file, header=None, names=["Chromosome", "End"])
    chromsizes.insert(1, "Start", 0)
    return chromsizes


def read_barcode_metrics(barcode_metrics_file: str, fragments_col: str) -> DataFrame:
    """
    Reads cellranger-arc output per_barcode_metrics.csv file
    barcode_metrics_file (str): a path to per_barcode_metrics.csv file
    fragments_col str: a name of the column with counts of fragments per barcode
    """
    barcode_metrics = read_csv(barcode_metrics_file, usecols=["barcode", fragments_col])
    return barcode_metrics


def validate_celltype_columns(
    celltype_columns: List[str], reference_columns: List[str]
) -> None:
    """
    Check if celltype dataframe contains all necessary columns
    celltype_df (DataFrame): a list of columns
    reference_columns (List[str]): a list of expected columns
    """
    if CELLTYPE_COLUMNS.difference(celltype_columns):
        raise InvalidColumnName(
            "Cell-type annotation file must have `sample_id`, `barcode`, `celltype` columns"
        )


def read_celltype_annotation(sample_id: str, celltype_file: str) -> DataFrame:
    """
    Read cell-type annotation file to pandas DataFrame
    sample_id (str): sample identifier
    chromsizes_file (str): path to the file with chromosome lengths from the UCSC databases
    """
    # read chromsizes file
    celltypes = read_csv(celltype_file)
    # check if DataFrame contains all necessary rows
    validate_celltype_columns(celltypes.columns, CELLTYPE_COLUMNS)
    # subsample sample_id rows only
    celltypes = celltypes[celltypes.sample_id == sample_id].copy()
    return celltypes


def filter_nan_entries(sample_id: str, celltypes: DataFrame, dropna: bool) -> None:
    """
    Checks if there are any NaN values in celltypes dataframe for sample_id
    sample_id (str): sample identifier
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    """
    nan_etries = celltypes.isna().values.any()
    if nan_etries:
        if dropna:
            logging.warning(
                f"Cell-type annotation file contains NaN values for sample_id={sample_id}"
            )
            logging.info(
                f"Removing NaN values from celltype annotation for {sample_id}"
            )
            celltypes.dropna(inplace=True)
        else:
            raise DataError(
                f"Cell-type annotation file contains NaN values for sample_id={sample_id}"
            )


def find_zero_fragment_celltypes(
    sample_id: str,
    celltypes: DataFrame,
    barcode_metrics: DataFrame,
    celltype_col: str,
    fragments_col: str,
) -> List[str]:
    """
    Find all celltypes for which there are no fragments found
    sample_id (str): sample identifier
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    barcode_metrics (str): a DataFrame with barcode metrics
    celltype_col (str): a name of the column with celltype labels in celltypes DataFrame
    fragments_col str: a name of the column with counts of fragments per barcode
    """
    # join two DataFrames
    joined_df = celltypes.merge(barcode_metrics, left_on="barcode", right_on="barcode")
    # calculate a total number of fragments for each celltype
    total_fragments = joined_df.groupby(celltype_col).sum(numeric_only=True)
    total_fragments.to_csv(f"{sample_id}_celltypes_fragment_counts.csv")
    # get celltypes with zero fragments
    zero_fragments_celltypes = total_fragments[
        total_fragments[fragments_col] == 0
    ].index.to_list()
    return zero_fragments_celltypes


def drop_empty_fragment(
    sample_id: str,
    celltypes: DataFrame,
    zero_fragments_celltypes: List[str],
    celltype_col: str,
):
    """
    Drops celltypes with zero fragments from celltype DataFrame
    sample_id (str): sample identifier
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    zero_fragments_celltypes (List[str]): a list of celltypes with zero fragments
    celltype_col (str): a name of the column with celltype labels in celltypes DataFrame
    """
    # log this step
    celltype_string = celltypes_string = ",".join(zero_fragments_celltypes)
    logging.info(
        f"Dropping the following celltypes from the {sample_id}'s annotation DataFrame as no fragments where found: {celltypes_string}"
    )
    # drop the celltypes
    zero_fragment_cells = celltypes[
        celltypes[celltype_col].isin(zero_fragments_celltypes)
    ].index
    celltypes.drop(zero_fragment_cells, inplace=True)


def filter_empty_fragments(
    sample_id: str,
    celltypes: DataFrame,
    barcode_metrics: DataFrame,
    celltype_col: str,
    fragments_col: str,
    skip_empty_fragments: bool,
) -> DataFrame:
    """
    Filter celltypes with zero fragments
    sample_id (str): sample identifier
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    barcode_metrics (str): a DataFrame with barcode metrics
    celltype_col (str): a name of the column with celltype labels in celltypes DataFrame
    fragments_col str: a name of the column with counts of fragments per barcode
    skip_empty_fragments: if true removes all celltypes with zero fragments from DataFrame
    """
    zero_fragments_celltypes = find_zero_fragment_celltypes(
        sample_id, celltypes, barcode_metrics, celltype_col, fragments_col
    )
    if zero_fragments_celltypes:
        celltypes_string = ",".join(zero_fragments_celltypes)
        if skip_empty_fragments:
            logging.warning(
                f"Sample {sample_id} containts zero fragments for the following celltypes: {celltypes_string}"
            )
            drop_empty_fragment(
                sample_id, celltypes, zero_fragments_celltypes, celltype_col
            )
        else:
            raise ValueError(
                f"Sample {sample_id} containts zero fragments for the following celltypes: {celltypes_string}"
            )


def main():
    # set up logger
    setup_logging()

    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # read files to DataFrame
    chromsizes = read_chromsizes(args.chromsizes)
    celltypes = read_celltype_annotation(args.sample_id, args.celltype_annotation)
    barcode_metrics = read_barcode_metrics(args.barcode_metrics, args.fragments_col)

    # filter celltypes DataFrame
    filter_nan_entries(args.sample_id, celltypes, args.dropna)
    filter_empty_fragments(
        args.sample_id,
        celltypes,
        barcode_metrics,
        args.celltype_col,
        args.fragments_col,
        args.skip_empty_fragments,
    )

    # get fragments file path dict
    fragments = {args.sample_id: args.fragments}

    # make output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # run pseudobulk generation
    bw_paths, bed_paths = export_pseudobulk(
        input_data=celltypes,
        variable=args.celltype_col,
        sample_id_col=args.sample_id_col,
        chromsizes=chromsizes,
        bed_path=args.output_dir,
        bigwig_path=args.output_dir,
        path_to_fragments=fragments,
        n_cpu=args.cpus,
        normalize_bigwig=True,
    )


if __name__ == "__main__":
    main()
