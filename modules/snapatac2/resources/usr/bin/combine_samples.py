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
        "--samples_file", 
        type=str, 
        metavar="<val>", 
        help="Path to csv file that contains three columns: sample_id, path h5ad, path to gene_matrix h5ad",
    )
    
    parser.add_argument(
        "--genome",
        metavar="<file>",
        type=str,
        help="Specify genome to use. One of hg38, hg37, mm10, mm39, GRCh37, GRCh38, GRCm38, GRCm39. GRCh38 is used by default.",
    )
    
    parser.add_argument(
        "--n_features",
        metavar="<val>",
        type=int,
        help="Number of features to use (for doublets)",
    )
    
    parser.add_argument(
        "--postprocess",
        action = 'store_true',
        default = False,
        help="Whether to do dim reduction, batch correction (along samples) and clustering"
    )

    return parser


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()
    
    genome = getattr(snap.genome,args.genome)
    
    samples = pd.read_csv(args.samples_file,header=None)
    
    data = snap.AnnDataSet(
      adatas = [(samples.iloc[i,0],samples.iloc[i,1]) for i in range(samples.shape[0])],
      filename='combined.h5ad')
      
      
    unique_cell_ids = [sa + ':' + bc for sa, bc in zip(data.obs['sample'], data.obs_names)]
    data.obs_names = unique_cell_ids
    snap.pp.select_features(data,n_features=args.n_features)
    
    # save and close
    gene_matrix = snap.pp.make_gene_matrix(data,genome)
    gene_matrix.write_h5ad("combined_gene_matrix.h5ad",compression='gzip')
    
    
    # optional dim reduction and batch correction
    if(args.postprocess):
      snap.tl.spectral(data)
      snap.tl.umap(data)
    
      snap.pp.mnc_correct(data, batch="sample")
      snap.pp.harmony(data, batch="sample", max_iter_harmony=20)
    
      snap.tl.umap(data, use_rep="X_spectral_mnn")
      snap.tl.umap(data, use_rep="X_spectral_harmony")
    
      snap.pp.knn(data, use_rep="X_spectral_harmony")
      snap.tl.leiden(data)

    data.close()
    
if __name__ == "__main__":
    main()
