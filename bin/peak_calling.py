#!/usr/bin/env python
import argparse
import os
from tempfile import TemporaryDirectory
from typing import Dict
from pycisTopic.pseudobulk_peak_calling import peak_calling


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Performs pseudobulk peak calling with MACS2. It requires to have MACS2 installed (https://github.com/macs3-project/MACS)"
    )
    parser.add_argument(
        "--macs_path",
        type=str,
        metavar="<file>",
        help="Path to MACS binary (e.g. /xxx/MACS/xxx/bin/macs2)",
        default="macs2"
    )
    parser.add_argument(
        "--bed_path",
        metavar="<dir>",
        type=str,
        help="Specify a path to the directory with pseudobulk celltype.fragments.tsv.gz files"
    )
    parser.add_argument(
        "--output_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to the output directory"
    )
    parser.add_argument(
        "--specie",
        metavar="<val>",
        type=str,
        help="Effective genome size which is defined as the genome size which can be sequenced. Possible values: 'hs', 'mm', 'ce' and 'dm'. Default: 'hs'",
        default="hs"
    )
    parser.add_argument(
        "--cpus", metavar="<num>", type=int, help="Specify a number of cpu cores to use"
    )
    parser.add_argument(
        "--input_format",
        metavar="<val>",
        type=str,
        help="Format of tag file can be ELAND, BED, ELANDMULTI, ELANDEXPORT, SAM, BAM, BOWTIE, BAMPE, or BEDPE. Default is AUTO which will allow MACS to decide the format automatically. Default: BEDPE",
        default="BEDPE"
    )
    parser.add_argument(
        "--shift",
        metavar="<val>",
        type=int,
        help="To set an arbitrary shift in bp. For finding enriched cutting sites (such as in ATAC-seq) a shift of 73 bp is recommended. Default: 73",
        default=73
    )
    parser.add_argument(
        "--extend_read_size",
        metavar="<val>",
        type=int,
        help="To extend reads in 5->3 direction to fix-sized fragment. For ATAC-seq data, a extension of 146 bp is recommended. Default: 146",
        default=146
    )
    parser.add_argument(
        "--keep_duplicates",
        metavar="<val>",
        type=str,
        help="Whether to keep duplicate tags at te exact same location. Default: all",
        default="all"
    )
    parser.add_argument(
        "--q_value_cutoff",
        metavar="<val>",
        type=float,
        help="The q-value (minimum FDR) cutoff to call significant regions. Default: 0.05",
        default=0.05
    )
    parser.add_argument(
        "--skip_empty_peaks",
        help="If specified skips celltypes with no peaks found",
        action='store_true'
    )

    return parser


def make_bedpath_dict(bed_path: str) -> Dict[str, str]:
    """
    The function assumes that bed files are saved in the form celltype.fragments.tsv.gz in the folder.
    It creates the Dict[celltype, path].
    bed_path (str): path to the directory where all bed files are stored
    """
    bedpath_dict = dict()
    file_list = os.listdir(bed_path)
    for filename in file_list:
        # remove suffix in filename
        celltype = filename.replace(".fragments.tsv.gz", "")
        # write to the dict
        bedpath_dict[celltype] = os.path.join(bed_path, filename)
    return bedpath_dict


def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # make a bedpath dict
    bedpath_dict = make_bedpath_dict(args.bed_path)

    # run peak calling on pseudobulks
    # WARN: Here I explicitly make a temporary directory in temp folder with specified name because
    # of ray.init error (Issue #7724) that forbids a path length more than 107 bytes
    # https://github.com/ray-project/ray/issues/7724
    # https://github.com/python/cpython/issues/93852
    with TemporaryDirectory(prefix="pk", dir="/tmp") as temp_dir:
        narrow_peak_dict = peak_calling(
            macs_path=args.macs_path,
            bed_paths=bedpath_dict,
            outdir=args.output_dir,
            genome_size=args.specie,
            n_cpu=args.cpus,
            input_format=args.input_format,
            shift=args.shift,
            ext_size=args.extend_read_size,
            keep_dup=args.keep_duplicates,
            q_value=args.q_value_cutoff,
            skip_empty_peaks=args.skip_empty_peaks,
            _temp_dir=temp_dir,
        )


if __name__ == "__main__":
    main()
