#!/usr/bin/env python3

import argparse
import snapatac2 as snap
import pandas as pd
import scanpy as sc
import numpy as np


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Creates snapATAC2 h5ad from fragments, calculates tsse, doublet score, filter cellby thrs and generates gene score matrix"
    )
    parser.add_argument(
        "--sample_id", 
        type=str, 
        metavar="<val>", 
        help="Sample identificator",
    )
    parser.add_argument(
        "--fragments",
        metavar="<file>",
        type=str,
        help="Specify a path to the fragments.tsv.gz",
    )
    
    parser.add_argument(
        "--genome",
        metavar="<file>",
        type=str,
        help="Specify genome to use. One of hg38, hg37, mm10, mm39, GRCh37, GRCh38, GRCm38, GRCm39. GRCh38 is used by default."
    )
    
    parser.add_argument(
        "--min_counts",
        metavar="<val>",
        type=int,
        help="Minimal numer of fragments per cell",
        default = 5000,
    )
    parser.add_argument(
        "--max_counts",
        metavar="<val>",
        type=int,
        help="Maximal numer of fragments per cell",
        default = 1000000,
    )
    parser.add_argument(
        "--min_tsse",
        metavar="<val>",
        type=float,
        help="Threshold for cell TSS enrichment score",
        default = 10,
    )
    
    parser.add_argument(
        "--n_features",
        metavar="<val>",
        type=int,
        help="Number of features to use (for doublets)",
    )
    return parser


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()
    
    genome = getattr(snap.genome,args.genome)
    
    data = snap.pp.import_fragments(
      args.fragments,
      chrom_sizes=genome,
      file=args.sample_id+".h5ad",
      sorted_by_barcode=False
    )
    
    snap.metrics.tsse(data,genome)
    snap.pp.filter_cells(data,
                         min_counts = args.min_counts,
                         min_tsse = args.min_tsse,
                         max_counts = args.max_counts)
                         
    if data.shape[0] == 0:
      raise Exception("There are no cells left in '"+args.sample_id+"', cannot proceed. Consider removing sample from a list.")
    
    snap.pp.add_tile_matrix(data)
    snap.pp.select_features(data,
                            n_features = args.n_features)
    snap.pp.scrublet(data)
    snap.pp.filter_doublets(data)
    
    #gene_matrix = snap.pp.make_gene_matrix(data,genome)
    #gene_matrix.write_h5ad(args.sample_id+"_gene_matrix.h5ad",compression='gzip')
    
    data.close()
    


if __name__ == "__main__":
    main()
