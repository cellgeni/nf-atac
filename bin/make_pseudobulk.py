#!/usr/bin/env python
import argparse
import os
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


def read_celltype_annotation(sample_id: str, celltype_file: str) -> DataFrame:
    """
    Read cell-type annotation file to pandas DataFrame
    chromsizes_file (str): path to the file with chromosome lengths from the UCSC databases
    """
    # read chromsizes file
    celltypes = read_csv(celltype_file)
    # check if DataFrame contains all necessary rows
    if CELLTYPE_COLUMNS.difference(celltypes.columns):
        raise InvalidColumnName(
            "Cell-type annotation file must have `sample_id`, `barcode`, `celltype` columns"
        )
    # subsample sample_id rows only
    celltypes = celltypes[celltypes.sample_id == sample_id].copy()
    # check if there are any NaN values
    if celltypes.isna().values.any():
        raise DataError(
            f"Cell-type annotation file contains NaN values for sample_id={sample_id}"
        )
    return celltypes


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # read files to DataFrame
    chromsizes = read_chromsizes(args.chromsizes)
    celltypes = read_celltype_annotation(args.sample_id, args.celltype_annotation)

    # get fragments file path dict
    fragments = {args.sample_id: args.fragments}

    # make output directory
    os.makedirs("output", exist_ok=True)

    # run pseudobulk generation
    bw_paths, bed_paths = export_pseudobulk(
        input_data=celltypes,
        variable="celltype",
        sample_id_col="sample_id",
        chromsizes=chromsizes,
        bed_path="output",
        bigwig_path="output",
        path_to_fragments=fragments,
        n_cpu=args.cpus,
        normalize_bigwig=True,
    )


if __name__ == "__main__":
    main()
