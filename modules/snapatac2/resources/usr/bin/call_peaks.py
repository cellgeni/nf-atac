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
        "--h5ad_file", 
        type=str,  
        help="Path to combined snapatac2 h5ad file",
    )
    parser.add_argument(
        "--celltype_file", 
        type=str,  
        help="Path to combined celltype annotation csv file with three columns: library_id,barcode,celltype",
    )
    parser.add_argument(
        "--genome",
        type=str,
        help="Specify genome to use. One of hg38, hg37, mm10, mm39, GRCh37, GRCh38, GRCm38, GRCm39. GRCh38 is used by default."
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        help="Number of processes"
    )
    return parser


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()
    
    genome = getattr(snap.genome,args.genome)
    
    data = snap.read_dataset(args.h5ad_file)
    celltypes = pd.read_csv(args.celltype_file)
    celltypes.index = celltypes.iloc[:,0]+":"+celltypes.iloc[:,1]
    celltypes.iloc[:,2] = celltypes.iloc[:,2].str.replace('/','_')

    cmn = list(set(celltypes.index).intersection(set(data.obs_names)))
    
    data,ord = data.subset(obs_indices=cmn,out='subset_adatas')
    data.obs['celltype'] = celltypes.loc[data.obs_names,:].iloc[:,2]

    snap.tl.macs3(data, groupby='celltype', replicate='sample',n_jobs=args.n_jobs)
    merged_peaks = snap.tl.merge_peaks(data.uns['macs3'], chrom_sizes=genome)
    
    peak_mat = snap.pp.make_peak_matrix(data, use_rep=merged_peaks['Peaks'])
    peak_mat.write_h5ad('subset_adatas/peak_mat.h5ad')
    merged_peaks.write_csv('subset_adatas/merged_peaks.csv')

    data.close()
    
if __name__ == "__main__":
    main()
