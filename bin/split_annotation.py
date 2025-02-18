#!/usr/bin/env python3

import os
import logging
import argparse
import collections
from typing import List, Dict
from pandas import read_csv, concat, DataFrame, Series
from colored_logger import setup_logging
from pandas.errors import InvalidColumnName, DataError


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Script validates sample and annotation tables and splits annotation table into separate celltypes"
    )
    parser.add_argument(
        "--sample_id",
        type=str,
        metavar="<val>",
        nargs="+",
        help="Specify sample ids in celltype annotation file you want to process",
    )
    parser.add_argument(
        "--celltype_annotation",
        metavar="<file>",
        type=str,
        help="Specify a path to the cell-type annotation file",
    )
    parser.add_argument(
        "--barcode_metrics",
        type=str,
        metavar="<file>",
        nargs="+",
        help="Specify Cellranger's per_barcode_metrics.csv files in the same order as sample ids",
    )
    parser.add_argument(
        "--sample_table",
        type=str,
        metavar="<file>",
        help="Specify sample table file",
    )
    parser.add_argument(
        "--output_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to the output directory",
        default="output",
    )
    parser.add_argument(
        "--updated_sample_table",
        metavar="<file>",
        type=str,
        help="Specify a name for the updated sample table",
        default="updated_sample_table.csv",
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
        help="Specify a name for sample_id column in annotation and sample table",
        default="sample_id",
    )
    parser.add_argument(
        "--barcode_col",
        metavar="<val>",
        type=str,
        help="Specify a name for barcode column in annotation file",
        default="barcode",
    )
    parser.add_argument(
        "--fragments_col",
        metavar="<val>",
        type=str,
        help="Specify a name for fragments column in 'per_barcode_metrics.csv' file",
        default="atac_fragments",
    )
    parser.add_argument(
        "--filedir_col",
        metavar="<val>",
        type=str,
        help="Specify a name for a cellranger-arc output directory column in sample table",
        default="filedir",
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
        "--logfile",
        metavar="<file>",
        type=str,
        help="Specify a log file name",
        default="splitcelltypes.log",
    )
    parser.add_argument(
        "--fragments_celltype_x_sample",
        metavar="<file>",
        type=str,
        help="Specify a name for the fragments_celltype_x_sample file",
        default="fragments_celltype_x_sample.csv",
    )
    parser.add_argument(
        "--fragments_per_celltype",
        metavar="<file>",
        type=str,
        help="Specify a name for the fragments_per_celltype file",
        default="fragments_per_celltype.csv",
    )
    return parser


def validate_sample_table_columns(
    sample_table_columns: List[str], sample_id_col: str, filedir_col: str
) -> None:
    """
    Check if sample dataframe contains all necessary columns
    Args:
        sample_table_columns (List[str]): a list of columns in sample table
        sample_id_col (str): a name for sample_id column in sample files
        filedir_col (str): a name for cellranger-arc output directory column in annotation and sample files

    Raises:
        InvalidColumnName: if column names differ from reference names
    """
    reference_columns = [sample_id_col, filedir_col]
    if reference_columns != sample_table_columns:
        raise InvalidColumnName(
            "Sample table file must have '{sample_id_col}' and '{filedir_col}' columns"
        )


def read_sample_table(filepath: str, sample_id_col: str, filedir_col: str) -> DataFrame:
    """
    Read and validate sample table
    Args:
        filepath (str): a path to sample table
        sample_id_col (str): a name for sample_id column in sample files
        filedir_col (str): a name for cellranger-arc output directory column in annotation and sample files

    Returns:
        DataFrame: sample table
    """
    # read to DataFrame
    sample_df = read_csv(filepath)
    # check if column names are valid
    df_columns = sample_df.columns.to_list()
    validate_sample_table_columns(df_columns, sample_id_col, filedir_col)
    return sample_df


def validate_sample_list(samples: List[str]) -> None:
    """
    Check if there are duplicated samples in the sample list

    Args:
        samples (List[str]): a list of sample ids
    Raises:
        DataError: if there is duplicates in the list of samples
    """
    if len(samples) != len(set(samples)):
        duplicated_samples = [
            sample
            for sample, count in collections.Counter(samples).items()
            if count > 1
        ]
        duplicated_samples_string = ",".join(duplicated_samples)
        raise DataError(
            f"Multiple entries of {duplicated_samples_string} in sample table"
        )


def validate_celltype_columns(
    celltype_columns: List[str], sample_id_col: str, celltype_col: str, barcode_col: str
) -> None:
    """
    Check if celltype DataFrame contains all necessary columns
    celltype_columns (List[str]):
    reference_columns (Set[str]): a list of expected columns

    Args:
        celltype_columns (List[str]): a list of columns in celltype annotation
        sample_id_col (str): a columns with sample id in annotation file
        celltype_col (str): a columns with celltype label in annotation file
        barcode_col (str): a columns with barcode sequence in annotation file

    Raises:
        InvalidColumnName: if any columns are missing
    """
    reference_columns = {sample_id_col, celltype_col, barcode_col}
    if reference_columns.difference(celltype_columns):
        raise InvalidColumnName(
            "Cell-type annotation file must have '{sample_id_col}', '{celltype_col}' and '{barcode_col}' columns"
        )


def filter_nan_entries(celltypes: DataFrame, dropna: bool) -> None:
    """
    Filters NaN values in celltypes DataFrame for sample_id
    Args:
        celltypes (DataFrame): a DataFrame with celltype annotation
        dropna (bool): if true drops NaN values from celltype annotation

    Raises:
        DataError: if dropna=False and annotation file contains NaN value
    """
    nan_entries = celltypes.isna().values.any()
    if nan_entries:
        if dropna:
            logging.warning(
                f"Cell-type annotation file contains NaN values. Removing them..."
            )
            celltypes.dropna(inplace=True)
        else:
            raise DataError(f"Cell-type annotation file contains NaN values")
    return celltypes


def subsample_celltypes(
    celltypes: DataFrame, samples: List[str], sample_id_col: str
) -> DataFrame:
    """
    Subsample samples from celltypes annotation
    Args:
        celltypes (DataFrame): a DataFrame with celltype annotation
        samples (List[str]): a list of samples from sample table
        sample_id_col (str): a name for sample_id column in annotation file

    Raises:
        DataError: if there are no such such samples in the annotation file

    Returns:
        DataFrame: subsampled celltype annotation DataFrame
    """
    # check if all samples are in the celltype annotation DataFrame
    samples_in_annotation = celltypes[sample_id_col].unique().tolist()
    difference = set(samples).difference(samples_in_annotation)
    if difference:
        sample_string = ",".join(difference)
        raise DataError(
            f"There are no {sample_string} samples in celltype annotation file"
        )
    # subsample entries from sample set
    celltypes_subsample = celltypes[celltypes[sample_id_col].isin(samples)].copy()
    return celltypes_subsample


def read_celltype_annotation(
    filepath: str,
    samples: List[str],
    sample_id_col: str,
    celltype_col: str,
    barcode_col: str,
    dropna: bool,
) -> DataFrame:
    """
    Read cell-type annotation file to pandas DataFrame
    Args:
        filepath (str): a path to celltype annotation file
        samples (List[str]): a list of samples from sample table
        sample_id_col (str): a name for sample_id column in annotation file
        celltype_col (str): a name for celltype column in annotation file
        barcode_col (str): a name for barcode column in annotation file
        dropna (bool): if true drops NaN values from celltype annotation

    Returns:
        DataFrame: cell-type annotation for samples in sample list
    """
    # read chromsizes file
    celltypes = read_csv(filepath)
    # check if DataFrame contains all necessary rows
    validate_celltype_columns(
        celltypes.columns, sample_id_col, celltype_col, barcode_col
    )
    # filter NaN entries
    celltypes_filtered = filter_nan_entries(celltypes, dropna)
    # leave only samples from sample_table
    celltypes_subsample = subsample_celltypes(
        celltypes_filtered, samples, sample_id_col
    )
    return celltypes_subsample


def read_barcode_metrics(
    barcode_metrics_file: str, barcode_col: str, fragments_col: str
) -> DataFrame:
    """
    Reads cellranger-arc output per_barcode_metrics.csv file
    Args:
        barcode_metrics_file (str): a path to per_barcode_metrics.csv file
        barcode_col str: a name of the column with barcode names
        fragments_col str: a name of the column with counts of fragments per barcode

    Returns:
        DataFrame: contains fragment number for each barcode
    """
    barcode_metrics = read_csv(
        barcode_metrics_file, usecols=[barcode_col, fragments_col]
    )
    return barcode_metrics


def get_sample_fragments_per_celltype(
    sample: str,
    barcode_metrics_df: DataFrame,
    celltypes: DataFrame,
    sample_id_col: str,
    celltype_col: str,
    barcode_col: str,
    fragments_col: str,
) -> DataFrame:
    """
    Generates a DataFrame of fragment number per celltype for sample
    Args:
        sample (str): _description_
        barcode_metrics_df (DataFrame): a DataFrame with fragment number for each barcode
        celltypes (DataFrame): a DataFrame with celltype annotation
        sample_id_col (str): a name for sample_id column in annotation file
        celltype_col (str): a name for celltype column in annotation file
        barcode_col (str): a name for barcode column in annotation and barcode metrics files
        fragments_col (str): a name for fragments number column in barcode metrics file

    Returns:
        DataFrame: contains a fragment number per celltype, shape=(n_celltypes, 1)
    """
    # subsample celltype DataFrame
    celltypes_sub = celltypes[celltypes[sample_id_col] == sample].copy()
    # join celltype DataFrames
    joined_df = celltypes_sub.merge(
        barcode_metrics_df, left_on=barcode_col, right_on=barcode_col
    )
    # calculate a total number of fragments per celltype
    total_fragments = (
        joined_df.groupby(celltype_col)
        .aggregate({fragments_col: "sum"})
        .rename({fragments_col: sample}, axis=1)
    )
    return total_fragments


def get_fragment_counts(
    sample_to_metrics: Dict[str, DataFrame],
    celltypes: DataFrame,
    sample_id_col: str,
    celltype_col: str,
    barcode_col: str,
    fragments_col: str,
) -> DataFrame:
    """

    Args:
        sample_to_metrics (Dict[str, DataFrame]): a Dict with barcode metrics DataFrames for each sample
        celltypes (DataFrame): a DataFrame with celltype annotation
        sample_id_col (str): a name for sample_id column in annotation file
        celltype_col (str): a name for celltype column in annotation file
        barcode_col (str): a name for barcode column in annotation and barcode metrics files
        fragments_col (str): a name for fragments number column in barcode metrics file

    Returns:
        DataFrame: a number of fragments per celltype for each sample, shape=(n_celltypes, n_samples)
    """
    # calculate total fragments per celltype for each sample
    sample_fragments_per_celltype_list = [
        get_sample_fragments_per_celltype(
            sample,
            barcode_metrics,
            celltypes,
            sample_id_col,
            celltype_col,
            barcode_col,
            fragments_col,
        )
        for sample, barcode_metrics in sample_to_metrics.items()
    ]
    # convert to DataFrame
    fragments_celltype_x_sample = (
        concat(sample_fragments_per_celltype_list, axis=1).fillna(0).astype("int64")
    )
    return fragments_celltype_x_sample


def get_fragments_per_sample(
    sample_to_metrics: Dict[str, DataFrame], fragments_col: str
) -> DataFrame:
    """
    Creates a DataFrame with total number of fragments for each sample
    Args:
        sample_to_metrics (Dict[str, DataFrame]): a Dict with barcode metrics DataFrames for each sample
        fragments_col (str): a name for fragments number column in barcode metrics file

    Returns:
        DataFrame: a total number of fragments for each sample, shape=(n_samples, 1)
    """
    sample_fragments_list = [
        metrics[fragments_col].sum() for sample, metrics in sample_to_metrics.items()
    ]
    fragments_per_sample_df = DataFrame(
        sample_fragments_list, columns=["fragments"], index=sample_to_metrics.keys()
    )
    return fragments_per_sample_df


def get_fragments_per_celltype(fragments_celltype_x_sample: DataFrame) -> Series:
    """
    Get fragments per celltype from fragments_celltype_x_sample DataFrame
    Args:
        fragments_celltype_x_sample (DataFrame): fragment counts table of shape=(n_celltypes, n_samples)

    Returns:
        Series: fragments per celltype
    """
    fragments_per_celltype = fragments_celltype_x_sample.sum(axis=1)
    return fragments_per_celltype


def drop_zero_fragment_celltypes(
    fragments_per_celltype: Series, celltypes: DataFrame, celltype_col: str
) -> DataFrame:
    """
    Drops celltypes with zero fragments from celltype DataFrame
    Args:
        fragments_per_celltype (Series): a Series with total number of fragments per celltype
        celltypes (DataFrame): a DataFrame with celltype annotation
        celltype_col (str): a name for celltype column in annotation file

    Returns:
        DataFrame: filtered celltype annotation
    """
    # find celltypes with zero fragments
    zero_fragment_celltypes = fragments_per_celltype[
        fragments_per_celltype == 0
    ].index.to_list()
    if zero_fragment_celltypes:
        celltype_string = ",".join(zero_fragment_celltypes)
        # log this step
        logging.warning(
            f"No fragments for {celltype_string} celltypes. Dropping them from celltype annotation file..."
        )
        # drop the celltypes
        zero_fragment_rows = celltypes[
            celltypes[celltype_col].isin(zero_fragment_celltypes)
        ].index
        celltypes = celltypes.drop(zero_fragment_rows).copy()
    return celltypes


def split_celltypes(celltypes_df: DataFrame, celltype_col: str, output_dir: str):
    """
    Split celltype annotation in separate files for each celltype
    Args:
        celltypes_df (DataFrame): a DataFrame with celltype annotation
        celltype_col (str): a name for celltype column in annotation file
        output_dir (str): an output directory name
    """
    unique_celltypes = celltypes_df[celltype_col].unique().tolist()
    for celltype in unique_celltypes:
        # subsample df for separate celltypes
        celltype_df = celltypes_df[celltypes_df[celltype_col] == celltype].copy()
        # save to csv
        celltype_filename = os.path.join(
            output_dir, f"{celltype.replace(' ', '_')}.csv"
        )
        celltype_df.to_csv(celltype_filename, index=False)


def update_sample_table(
    sample_table: DataFrame,
    sample_to_metrics: Dict[str, DataFrame],
    sample_id_col: str,
    fragments_col: str,
) -> DataFrame:
    """
    Add fragment count to sample table
    Args:
        sample_table (DataFrame): a sample table
        sample_to_metrics (Dict[str, DataFrame]): a Dict with barcode metrics DataFrames for each sample
        sample_id_col (str): a name for sample_id column in sample table
        fragments_col (str): a name for fragments number column in barcode metrics file

    Returns:
        DataFrame: Updated sample table with fragment counts
    """
    # set sample_id as index for sample_table
    updated_sample_table = sample_table
    # calculate fragments per sample
    fragments_per_sample = get_fragments_per_sample(
        sample_to_metrics, fragments_col
    ).reset_index(names=sample_id_col)
    # # join sample table and fragments per sample
    updated_sample_table = sample_table.merge(
        fragments_per_sample, on=sample_id_col
    ).set_index(sample_id_col)
    return updated_sample_table


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # set up logger
    setup_logging(args.logfile)

    # check if there are duplicates in sample list
    validate_sample_list(args.sample_id)

    # read sample table
    sample_table = read_sample_table(
        args.sample_table, args.sample_id_col, args.filedir_col
    )

    # read celltype annotation
    celltypes = read_celltype_annotation(
        filepath=args.celltype_annotation,
        samples=args.sample_id,
        sample_id_col=args.sample_id_col,
        celltype_col=args.celltype_col,
        barcode_col=args.barcode_col,
        dropna=True,
    )

    # make a dictionary with barcode metrics DataFrames
    sample_to_metrics_path = dict(zip(args.sample_id, args.barcode_metrics))
    sample_to_metrics = {
        sample: read_barcode_metrics(metrics_path, args.barcode_col, args.fragments_col)
        for sample, metrics_path in sample_to_metrics_path.items()
    }

    # calculate fragment numbers
    fragments_celltype_x_sample = get_fragment_counts(
        sample_to_metrics,
        celltypes,
        args.sample_id_col,
        args.celltype_col,
        args.barcode_col,
        args.fragments_col,
    )
    fragments_per_celltype = get_fragments_per_celltype(fragments_celltype_x_sample)

    # update sample table
    updated_sample_table = update_sample_table(
        sample_table, sample_to_metrics, args.sample_id_col, args.fragments_col
    )

    # fillter celltypes with zero fragments and split celltype annotation in separate files
    os.mkdir(args.output_dir)
    celltypes = drop_zero_fragment_celltypes(
        fragments_per_celltype, celltypes, args.celltype_col
    )
    split_celltypes(celltypes, args.celltype_col, args.output_dir)

    # save fragment counts
    fragments_celltype_x_sample.to_csv(args.fragments_celltype_x_sample)
    fragments_per_celltype.index = fragments_per_celltype.index.str.replace(" ", "_")
    fragments_per_celltype.to_csv(args.fragments_per_celltype, header=False)
    updated_sample_table.to_csv(args.updated_sample_table)


if __name__ == "__main__":
    main()
