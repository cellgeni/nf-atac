#!/usr/bin/env python3

import os
import logging
import argparse
from typing import List, Set
from pandas import read_csv, DataFrame
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
        "--sample_table",
        metavar="<file>",
        type=str,
        help="Specify a path to the sample file",
    )
    parser.add_argument(
        "--celltype_annotation",
        metavar="<file>",
        type=str,
        help="Specify a path to the cell-type annotation file",
    )
    parser.add_argument(
        "--output_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to the output directory",
        default="output",
    )
    parser.add_argument(
        "--filtered_sample_table",
        metavar="<file>",
        type=str,
        help="Specify a name for the filtered sample table",
        default="filtered_sample_table.csv",
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
        help="Specify a name for sample_id column in annotation and sample files",
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
        "--filedir_col",
        metavar="<val>",
        type=str,
        help="Specify a name for a cellranger-arc output directory column in sample file",
        default="filedir",
    )
    parser.add_argument(
        "--dropna",
        help="If specified drops NaN values from celltype annotation",
        action="store_true",
    )
    parser.add_argument(
        "--drop_duplicated_samples",
        help="If specified drops duplicated samples from sample table",
        action="store_true",
    )
    parser.add_argument(
        "--logfile",
        metavar="<file>",
        type=str,
        help="Specify a log file name",
        default="splitcelltypes.log",
    )
    return parser


def validate_sample_table_columns(
    sample_table_columns: List[str], sample_id_col: str, filedir_col: str
) -> None:
    """
    Check if sample dataframe contains all necessary columns
    sample_table_columns (List[str]): a list of columns in sample table
    sample_id_col (str): a name for sample_id column in sample files
    filedir_col (str): a name for cellranger-arc output directory column in annotation and sample files
    """
    reference_columns = [sample_id_col, filedir_col]
    if reference_columns != sample_table_columns:
        raise InvalidColumnName(
            "Sample table file must have '{sample_id_col}' and '{filedir_col}' columns"
        )


def drops_duplicated_sample_ids(
    sample_df: DataFrame, sample_id_col: str, drop_duplicated_samples: bool
) -> DataFrame:
    """
    Remove dupicated sample ids
    sample_df (DataFrame): a DataFrame with sample table
    sample_id_col (str): a name for sample_id column in asample file
    drop_duplicated_samples (bool): if true drops duplicated samples from sample table
    """
    # find duplicated samples
    duplicated = sample_df[sample_id_col].duplicated(keep="first")
    if duplicated.any() and drop_duplicated_samples:
        # raise warning
        dup_samples_string = ",".join(
            sample_df[duplicated]["sample_id_col"].unique().tolist()
        )
        logging.warning(
            f"Multiple entries of {dup_samples_string} in sample table. Dropping duplicated..."
        )
        # drop duplicated samples
        sample_df_filtered = sample_df[~duplicated].copy()
        return sample_df_filtered
    elif duplicated.any() and not drop_duplicated_samples:
        # raise an error
        dup_samples_string = ",".join(
            sample_df[duplicated]["sample_id_col"].unique().tolist()
        )
        raise DataError(f"Multiple entries of {dup_samples_string} in sample table")
    return sample_df


def read_sample_table(
    filepath: str, sample_id_col: str, filedir_col: str, drop_duplicated_samples: bool
) -> DataFrame:
    """
    Read and validate sample table
    filepath (str): a path to sample table
    sample_id_col (str): a name for sample_id column in sample files
    filedir_col (str): a name for cellranger-arc output directory column in annotation and sample files
    drop_duplicated_samples (bool): if true drops duplicated samples from sample table
    """
    # read to DataFrame
    sample_df = read_csv(filepath)
    # check if column names are valid
    df_columns = sample_df.columns.to_list()
    validate_sample_table_columns(df_columns, sample_id_col, filedir_col)
    # remove duplicated sample ids
    sample_df = drops_duplicated_sample_ids(
        sample_df, sample_id_col, drop_duplicated_samples
    )
    return sample_df


def validate_celltype_columns(
    celltype_columns: List[str], sample_id_col: str, celltype_col: str, barcode_col: str
) -> None:
    """
    Check if celltype dataframe contains all necessary columns
    celltype_columns (List[str]): a list of columns in celltype annotation
    reference_columns (Set[str]): a list of expected columns
    """
    reference_columns = {sample_id_col, celltype_col, barcode_col}
    if reference_columns.difference(celltype_columns):
        raise InvalidColumnName(
            "Cell-type annotation file must have '{sample_id_col}', '{celltype_col}' and '{barcode_col}' columns"
        )


def filter_nan_entries(celltypes: DataFrame, dropna: bool) -> None:
    """
    Filters NaN values in celltypes dataframe for sample_id
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    dropna (bool): if true drops NaN values from celltype annotation
    """
    nan_etries = celltypes.isna().values.any()
    if nan_etries:
        if dropna:
            logging.warning(
                f"Cell-type annotation file contains NaN values. Removing them..."
            )
            celltypes.dropna(inplace=True)
        else:
            raise DataError(f"Cell-type annotation file contains NaN values")


def subsample_celltypes(
    celltypes: DataFrame, samples: Set[str], sample_id_col: str
) -> DataFrame:
    """
    Subsample samples from celltypes annotation
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    samples (Set[str]): a set of samples from sample table
    sample_id_col (str): a name for sample_id column in annotation file
    """
    # check if all samples are in the celltype annotation DataFrame
    samples_in_annotation = celltypes[sample_id_col].unique().tolist()
    difference = samples.difference(samples_in_annotation)
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
    samples: Set[str],
    sample_id_col: str,
    celltype_col: str,
    barcode_col: str,
    dropna: bool,
) -> DataFrame:
    """
    Read cell-type annotation file to pandas DataFrame
    filepath (str): a path to celltype annotation file
    samples (Set[str]): a set of samples from sample table
    sample_id_col (str): a name for sample_id column in annotation file
    celltype_col (str): a name for celltype column in annotation file
    barcode_col (str): a name for barcode column in annotation file
    dropna (bool): if true drops NaN values from celltype annotation
    """
    # read chromsizes file
    celltypes = read_csv(filepath)
    # check if DataFrame contains all necessary rows
    validate_celltype_columns(
        celltypes.columns, sample_id_col, celltype_col, barcode_col
    )
    # filter NaN entries
    filter_nan_entries(celltypes, dropna)
    # leave only samples from sample_table
    celltypes_subsample = subsample_celltypes(celltypes, samples, sample_id_col)
    return celltypes_subsample


def split_celltypes(celltypes_df: DataFrame, celltype_col: str, output_dir: str):
    """
    Split celltype annotation in separate files for each celltype
    celltypes_df (DataFrame): a dataframe with barcode/celltype annotation for sample_id
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


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # set up logger
    setup_logging(args.logfile)

    # read sample table and get sample list
    sample_table = read_sample_table(
        args.sample_table,
        args.sample_id_col,
        args.filedir_col,
        args.drop_duplicated_samples,
    )
    samples = set(sample_table[args.sample_id_col])

    # read celltype annotation
    celltypes = read_celltype_annotation(
        args.celltype_annotation,
        samples,
        args.sample_id_col,
        args.celltype_col,
        args.barcode_col,
        args.dropna,
    )

    # split celltype annotation in separate files
    os.mkdir(args.output_dir)
    split_celltypes(celltypes, args.celltype_col, args.output_dir)

    # save filtered sample table
    sample_table.to_csv(args.filtered_sample_table, index=False)


if __name__ == "__main__":
    main()
