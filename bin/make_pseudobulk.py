#!/usr/bin/env python
import argparse
import os
import sys
import logging
from typing import List, Dict
from colored_logger import setup_logging
from pandas import read_csv, read_table, DataFrame, concat
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
        "--chromsizes",
        metavar="<file>",
        type=str,
        help="Specify a path to the file with chromosome lengths from the UCSC databases",
    )
    parser.add_argument(
        "--barcode_metrics",
        metavar="<file>",
        nargs="+",
        help="Specify Cellranger's per_barcode_metrics.csv files in the same order as sample ids",
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
        "--fragments_col",
        metavar="<val>",
        type=str,
        help="Specify a name for fragments column in 'per_barcode_metrics.csv' file",
        default="atac_fragments",
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


def filter_empty_fragments(
    celltypes: DataFrame,
    barcode_metrics_list: List[str],
    celltype_col: str,
    fragments_col: str,
    skip_empty_fragments: bool,
) -> DataFrame:
    """
    Filter celltypes with zero fragments for all samples in celltypes
    celltypes (DataFrame): a dataframe with barcode/celltype annotation
    barcode_metrics (List[str]): a List of barcode metrics file paths
    celltype_col (str): a name of the column with celltype labels in celltypes DataFrame
    fragments_col str: a name of the column with counts of fragments per barcode
    skip_empty_fragments: if true removes all celltypes with zero fragments from DataFrame
    """
    celltype_filtered = list()
    for sample_id, barcode_metrics_path in barcode_metrics_list.items():
        # read barcode metrics file
        barcode_metrics = read_barcode_metrics(barcode_metrics_path, fragments_col)
        # subsample the celltype annotation file
        celltype_sample = celltypes[celltypes.sample_id == sample_id].copy()
        # filter annotation from celltypes with no fragments
        celltype_filtered_sample = filter_empty_fragments_for_sample(
            sample_id=sample_id,
            celltypes=celltype_sample,
            barcode_metrics=barcode_metrics,
            celltype_col=celltype_col,
            fragments_col=fragments_col,
            skip_empty_fragments=skip_empty_fragments,
        )
        celltype_filtered.append(celltype_filtered_sample)
    # concat filtered dataframes
    celltypes = concat(celltype_filtered, axis=0)
    return celltypes


def check_missing_samples(
    celltypes: DataFrame, sample_list: List[str], sample_id_col: str, celltype_name: str
):
    """
    Check if there are samples with no barcodes available
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


def fragments_sum(
    sample_id: str,
    celltypes: DataFrame,
    barcode_metrics: DataFrame,
    celltype_col: str,
    fragments_col: str,
    barcode_col: str,
    celltype_name: str,
) -> int:
    """
    Calculates a number of fragments for the sample
    sample_id (str): sample identifier
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    barcode_metrics (str): a DataFrame with barcode metrics
    celltype_col (str): a name of the column with celltype labels in celltypes DataFrame
    fragments_col (str): a name of the column with counts of fragments per barcode
    barcode_col (str): a name of the column with barcode names in celltype annotation and barcode metrics files
    celltype_name (str): a name of the celltype in annotation file
    """
    # join two DataFrames
    joined_df = celltypes.merge(
        barcode_metrics, left_on=barcode_col, right_on=barcode_col
    )
    # calculate a total number of fragments for celltype
    total_fragments = joined_df.groupby(celltype_col).sum(numeric_only=True)
    # get celltypes fragments num
    fragments_num = total_fragments.at[celltype_name, fragments_col]
    logging.info(f"There are {fragments_num} fragments found for {sample_id}")
    return fragments_num


def filter_empty_fragments_for_sample(
    sample_id: str,
    celltypes: DataFrame,
    barcode_metrics: DataFrame,
    celltype_col: str,
    fragments_col: str,
    barcode_col: str,
    celltype_name: str,
    skip_empty_fragments: bool,
) -> bool:
    """
    Whether we should drop the sample
    sample_id (str): sample identifier
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    barcode_metrics (str): a DataFrame with barcode metrics
    celltype_col (str): a name of the column with celltype labels in celltypes DataFrame
    fragments_col (str): a name of the column with counts of fragments per barcode
    barcode_col (str): a name of the column with barcode names in celltype annotation and barcode metrics files
    celltype_name (str): a name of the celltype in annotation file
    skip_empty_fragments (bool): if true drops fragments with zero fragments found
    """
    fragments_num = fragments_sum(
        sample_id,
        celltypes,
        barcode_metrics,
        celltype_col,
        fragments_col,
        barcode_col,
        celltype_name,
    )
    if fragments_num == 0 and skip_empty_fragments:
        logging.warning(
            f"Sample {sample_id} containts zero fragments for the '{celltype_name}' celltype. Dropping sample..."
        )
        return True
    elif fragments_num == 0 and not skip_empty_fragments:
        raise ValueError(
            f"Sample {sample_id} containts zero fragments for the '{celltype_name}' celltype"
        )
    return False


def drop_empty_samples(
    sample_list: List[str], celltypes: DataFrame, sample_id_col: str, celltype_name: str
) -> DataFrame:
    """
    Drops sample with zero fragments from celltype DataFrame
    sample_list (List[str]): a list of samples that we want to drop
    celltypes (DataFrame): a dataframe with barcode/celltype annotation for sample_id
    sample_id_col (str): a name of the column with sample_id labels in celltypes DataFrame
    celltype_name (str): a name of the celltype in annotation file
    """
    if sample_list:
        sample_string = ",".join(sample_list)
        # log this step
        logging.info(
            f"Dropping {sample_string} from the from {celltype_name} annotation file"
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
    barcode_metrics_dict: Dict[str, str],
    sample_id_col: str,
    celltype_col: str,
    fragments_col: str,
    barcode_col: str,
    celltype_name: str,
    skip_empty_fragments: bool,
) -> DataFrame:
    """
    Filter samples with zero fragments found
    celltypes (DataFrame): a dataframe with barcode/celltype annotation
    barcode_metrics_dict (Dict[str, str]): a Dict with paths to barcode metrics files
    sample_id_col (str): a name of the column with sample_id labels in celltypes DataFrame
    celltype_col (str): a name of the column with celltype labels in celltypes DataFrame
    fragments_col (str): a name of the column with counts of fragments per barcode
    barcode_col (str): a name of the column with barcode names in celltype annotation and barcode metrics files
    celltype_name (str): a name of the celltype in annotation file
    skip_empty_fragments (bool): if true drops fragments with zero fragments found
    """
    samples_to_drop = list()
    # get samples from annotation
    samples_in_annotation = celltypes[sample_id_col].unique().tolist()
    for sample_id in samples_in_annotation:
        # subsample celltypes
        celltype_subsample = celltypes[celltypes[sample_id_col] == sample_id].copy()
        # read barcode metrics
        barcode_metrics = read_barcode_metrics(
            barcode_metrics_dict[sample_id], fragments_col
        )
        # check wether we want to drop sample
        drop_sample = filter_empty_fragments_for_sample(
            sample_id,
            celltype_subsample,
            barcode_metrics,
            celltype_col,
            fragments_col,
            barcode_col,
            celltype_name,
            skip_empty_fragments,
        )
        # save sample_id to list
        if drop_sample:
            samples_to_drop.append(sample_id)
    # drop samples with zero fragments found
    filtered_celltypes = drop_empty_samples(
        samples_to_drop, celltypes, sample_id_col, celltype_name
    )
    return filtered_celltypes


def read_celltype_annotation(
    filepath: str,
    sample_list: List[str],
    barcode_metrics_dict: Dict[str, str],
    sample_id_col: str,
    celltype_col: str,
    fragments_col: str,
    barcode_col: str,
    skip_empty_fragments: bool,
) -> DataFrame:
    """
    Read cell-type annotation file to pandas DataFrame
    filepath (str): a path to celltype annotation
    sample_list (List[str]): a list of samples mentioned in sample table
    barcode_metrics_dict (Dict[str, str]): a Dict with paths to barcode metrics files
    sample_id_col (str): a name of the column with sample_id labels in celltypes DataFrame
    celltype_col (str): a name of the column with celltype labels in celltypes DataFrame
    fragments_col (str): a name of the column with counts of fragments per barcode
    barcode_col (str): a name of the column with barcode names in celltype annotation and barcode metrics files
    skip_empty_fragments (bool): if true drops fragments with zero fragments found
    """
    # read chromsizes file
    celltype_name = filepath.replace(".csv", "").replace("_", " ")
    celltypes = read_csv(filepath)
    # check if all samples are present in annotation
    check_missing_samples(celltypes, sample_list, sample_id_col, celltype_name)
    # check if there are samples with zero fragments for this celltype
    filtered_celltypes = filter_empty_fragments(
        celltypes,
        barcode_metrics_dict,
        sample_id_col,
        celltype_col,
        fragments_col,
        barcode_col,
        celltype_name,
        skip_empty_fragments,
    )
    # check if DataFrame is empty
    if filtered_celltypes.empty:
        logging.info(f"There is no fragments found for {celltype_name} celltype")
        raise SystemExit
    return filtered_celltypes


def main():
    # set up logger
    setup_logging()

    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # make a dictionaries for fragments and barcode metrics files
    fragments_dict = dict(zip(args.sample_id, args.fragments))
    barcode_metrics_dict = dict(zip(args.sample_id, args.barcode_metrics))

    # read files to DataFrame
    chromsizes = read_chromsizes(args.chromsizes)
    celltypes = read_celltype_annotation(
        args.celltype_annotation,
        args.sample_id,
        barcode_metrics_dict,
        args.sample_id_col,
        args.celltype_col,
        args.fragments_col,
        args.barcode_col,
        args.skip_empty_fragments,
    )

    # subsample samples present in celltype DataFrame
    filtered_sample_list = celltypes[args.sample_id_col].unique().tolist()
    filtered_fragments_dict = {
        sample_id: fragments_dict[sample_id] for sample_id in filtered_sample_list
    }

    # make output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.bigwig_dir, exist_ok=True)

    # run pseudobulk generation
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
    )


if __name__ == "__main__":
    main()
