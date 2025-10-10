#!/usr/bin/env python3

import sys
import logging
import argparse
import mudata as md
import anndata as ad


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
        description="Entangle multiome AnnData GEX and ATAC objects"
    )
    parser.add_argument(
        "--gex",
        metavar="<file>",
        type=str,
        help="Specify a path to the .h5ad files for Gene Expression data",
    )
    parser.add_argument(
        "--atac",
        metavar="<file>",
        type=str,
        help="Specify a path to the .h5ad files for ATAC data",
    )
    parser.add_argument(
        "--prefix",
        metavar="<file>",
        type=str,
        help="Specify a prefix for the output files",
    )
    parser.add_argument(
        "--logfile",
        metavar="<file>",
        type=str,
        default=".couple_multiome.log",
        help="Specify a path to the log file",
    )
    return parser


def main():
    # Parse arguments
    parser = init_parser()
    args = parser.parse_args()

    # Setup logging
    setup_logging(args.logfile)

    # Load AnnData objects
    logging.info("Loading Gene Expression .h5ad file")
    gex_adata = ad.read_h5ad(args.gex)
    logging.info(gex_adata)
    logging.info("Loading ATAC .h5ad file")
    atac_adata = ad.read_h5ad(args.atac)
    logging.info(atac_adata)

    # Subset barcodes to intersection
    logging.info("Subsetting barcodes to intersection")
    common_barcodes = gex_adata.obs_names.intersection(atac_adata.obs_names)
    gex_adata = gex_adata[common_barcodes].copy()
    atac_adata = atac_adata[common_barcodes].copy()
    logging.info("Number of common barcodes: %s", common_barcodes.shape[0])

    # Create MuData object
    logging.info("Creating MuData object")
    mdata = md.MuData({"gex": gex_adata, "atac": atac_adata}, obs_names=common_barcodes)
    logging.info(mdata)

    # Save all objects
    logging.info("Saving objects")
    gex_adata.write_h5ad(f"{args.prefix}_gex.h5ad")
    atac_adata.write_h5ad(f"{args.prefix}_atac.h5ad")
    mdata.write_h5mu(f"{args.prefix}_multiome.h5mu")


if __name__ == "__main__":
    main()
