#!/usr/bin/env python3

import argparse
import pandas as pd
import scanpy as sc


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Combines mts, obs and var csv into adata"
    )
    parser.add_argument(
        "--mtx", 
        type=str, 
        help="Path to mtx",
    )
    
    parser.add_argument(
        "--obs", 
        type=str, 
        help="Path to obs csv",
    )
    
    parser.add_argument(
        "--var",
        type=str,
        help="Path to var csv",
    )
    
    parser.add_argument(
        "--out_h5ad",
        type=str,
        help="name for output file",
    )

    return parser


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    adata = sc.read_mtx(args.mtx)
    adata.obs = pd.read_csv(args.obs, index_col=0)
    adata.var = pd.read_csv(args.var, index_col=0)

    adata.write_h5ad(args.out_h5ad)
    
    


        
if __name__ == "__main__":
    main()
