#!/usr/bin/env python3

import sys
import logging
import argparse
import anndata as an
import scanpy as sc
import pandas as pd


class ColoredFormatter(logging.Formatter):
    """
    Custom formatter to add colors to log messages
    """

    blue = "\n\033[94m"
    yellow = "\033[93m"
    red = "\033[91m"
    reset = "\033[0m"
    format_str = "%(levelname)s: %(message)s"

    FORMATS = {
        logging.INFO: blue + format_str + reset,
        logging.WARNING: yellow + format_str + reset,
        logging.ERROR: red + format_str + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def setup_logging(logfile=".python.log") -> None:
    """
    Setup logging configuration of the script
    """
    # a basic config to save logs to metadata.log
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        filename=logfile,
        filemode="w",
    )

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.WARNING)
    # tell the handler to use colored format
    console.setFormatter(ColoredFormatter())
    # add the handler to the root logger
    logging.getLogger("").addHandler(console)


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Script validates sample and annotation tables and splits annotation table into separate celltypes"
    )
    parser.add_argument(
        "--h5ad_file",
        type=str,
        metavar="<file>",
        help="Specify a path to .h5ad file",
        required=True,
    )
    parser.add_argument(
        "--sample_id",
        type=str,
        metavar="<str>",
        default=None,
        help="Specify sample_id in annotation file",
    )
    parser.add_argument(
        "--barcode_column",
        type=str,
        metavar="<str>",
        default="obs_names",
        help='Specify barcode column name in AnnData object to match with metadata. Specify "obs_names" if you wish to match with AnnData object obs_names. Default: "obs_names"',
    )
    parser.add_argument(
        "--metadata",
        type=str,
        metavar="<file>",
        help="Specify a path to celltype annotation .csv file",
        default=None,
    )
    parser.add_argument(
        "--how",
        type=str,
        metavar="<str>",
        choices=["left", "right", "outer", "inner"],
        default="inner",
        help='Type of merge to be performed. Options are "left", "right", "outer", "inner". Default: "inner"',
    )
    parser.add_argument(
        "--output",
        metavar="<path>",
        type=str,
        help="Specify a path to output .h5ad file",
    )
    parser.add_argument(
        "--logfile",
        metavar="<file>",
        type=str,
        default=".attach_celltypes.log",
        help="Specify a path to the log file",
    )

    return parser


def format_obs(adata: an.AnnData, sample_id: str) -> pd.DataFrame:
    """
    Format AnnData's .obs
    Args:
        adata (an.AnnData): AnnData object to format
        barcode_column (str): barcode column name in AnnData object (use "obs_names" to match with AnnData object obs_names)
        sample_id (str): sample_id to match with metadata
    Returns:
        str: sample_id column name in obs
        str: obs index name
        pd.DataFrame: formatted obs DataFrame
    """
    # Copy DataFrame to avoid modifying the original
    obs = adata.obs.copy()

    # Add sample_id to obs
    sample_id_column = (
        "sample_id" if sample_id not in obs.columns else "merge_sample_id"
    )
    obs[sample_id_column] = sample_id

    # Reset index
    obs_index_name = (
        obs.index.name
        if obs.index.name and obs.index.name not in obs.columns
        else "obs_names"
    )
    obs.reset_index(inplace=True, names=obs_index_name)
    return sample_id_column, obs_index_name, obs


def vaildate_obs(obs: pd.DataFrame, barcode_column: str, obs_index_name: str) -> None:
    """
    Validate obs DataFrame
    Args:
        obs (pd.DataFrame): obs DataFrame to validate
        barcode_column (str): barcode column name in AnnData object
        obs_index_name (str): obs index name
    Raises:
        ValueError: if barcode column contains NaN values
    Returns:
        str: final barcode column name
    """
    # Get final barcode column name
    barcode_column = obs_index_name if barcode_column == "obs_names" else barcode_column

    # Check if barcode column contains any NaN values
    if obs[barcode_column].isnull().any():
        error_msg = f"Barcode column '{barcode_column}' contains NaN values. Please check the input .h5 file"
        logging.error(error_msg)
        raise ValueError(error_msg)
    return barcode_column


def validate_metadata(metadata: pd.DataFrame, sample_id: str) -> None:
    """
    Validate metadata DataFrame
    Args:
        metadata (pd.DataFrame): metadata DataFrame to validate
        sample_id (str): sample_id to match with metadata
    Raises:
        ValueError: if any of the required columns are missing or contain NaN values
    """
    # Check if sample_id and barcode columns exist in metadata
    if ("barcode" not in metadata.columns) or ("sample_id" not in metadata.columns):
        error_msg = "Metadata file must contain 'barcode' and 'sample_id' columns. Please check the input metadata file"
        logging.error(error_msg)
        raise ValueError(error_msg)

    # Check if sample_id column contains sample_id given in the argument
    if sample_id not in metadata["sample_id"].unique():
        error_msg = f"Sample ID '{sample_id}' not found in metadata. Please check the input metadata file"
        logging.error(error_msg)
        raise ValueError(error_msg)

    # Check if barcode column contains any NaN values
    if metadata["barcode"].isnull().any():
        error_msg = (
            "Barcode column contains NaN values. Please check the input metadata file"
        )
        logging.error(error_msg)
        raise ValueError(error_msg)

    # Check if sample_id column contains any NaN values
    if metadata["sample_id"].isnull().any():
        error_msg = (
            "Sample ID column contains NaN values. Please check the input metadata file"
        )
        logging.error(error_msg)
        raise ValueError(error_msg)


def merge_metadata(
    obs: pd.DataFrame,
    metadata: pd.DataFrame,
    barcode_column: str,
    sample_id_column: str,
    sample_id: str,
    how: str = "inner",
) -> pd.DataFrame:
    """
    Merge metadata with obs DataFrame
    Args:
        obs (pd.DataFrame): obs DataFrame to merge
        metadata (pd.DataFrame): metadata DataFrame to merge
        barcode_column (str): barcode column name in obs
        sample_id_column (str): sample_id column name in obs
        sample_id (str): sample_id to match with metadata
    Returns:
        pd.DataFrame: merged DataFrame
    """
    # Merge metadata with obs
    obs_merged = obs.merge(
        metadata,
        left_on=[barcode_column, sample_id_column],
        right_on=["barcode", "sample_id"],
        how=how,
    )

    # Check if number of cells before and after merging is the same
    if obs_merged.shape[0] != obs.shape[0]:
        logging.warning(
            f"{sample_id}: Number of cells before merging: {obs.shape[0]}, Number of cells after merging: {obs_merged.shape[0]}"
        )
    return obs_merged


def main():
    # Parse arguments
    parser = init_parser()
    args = parser.parse_args()

    # Setup logging
    setup_logging(args.logfile)

    # Load AnnData object
    logging.info("Loading .h5 file to AnnData object")
    adata = sc.read_h5ad(args.h5ad_file)
    logging.info(adata)

    # Load cell metadata
    logging.info("Loading cell metadata")
    metadata = pd.read_csv(args.metadata)

    # Format AnnData's .obs
    sample_id_column, obs_index_name, obs = format_obs(adata, args.sample_id)

    # Validate obs and metadata
    logging.info("Validating metadata")
    barcode_column = vaildate_obs(obs, args.barcode_column, obs_index_name)
    validate_metadata(metadata, args.sample_id)

    # Merge metadata with obs
    logging.info("Merging metadata with obs")
    obs_merged = merge_metadata(
        obs, metadata, barcode_column, sample_id_column, args.sample_id, args.how
    )

    # Update AnnData's .obs with merged DataFrame
    logging.info("Updating AnnData's .obs with merged DataFrame")
    obs_merged.set_index(obs_index_name, inplace=True)
    adata_metadata = adata[obs_merged.index]
    adata_metadata.obs = obs_merged
    logging.info(adata_metadata)

    # Write AnnData object to .h5ad file
    logging.info("Writing AnnData object to .h5ad file")
    adata_metadata.write_h5ad(args.output)
    logging.info("Done!")


if __name__ == "__main__":
    main()
