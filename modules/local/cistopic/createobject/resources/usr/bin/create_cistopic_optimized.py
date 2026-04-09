#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import json
import logging
import os
import pickle
import sys
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import pyranges as pr
from pandas import DataFrame, read_parquet
from scipy import sparse

from pycisTopic.cistopic_class import CistopicObject
from pycisTopic.fragments import read_fragments_to_pyranges
from pycisTopic.qc import get_barcodes_passing_qc_for_sample


LOGGER_NAME = "cisTopic"


def get_logger() -> logging.Logger:
    logger = logging.getLogger(LOGGER_NAME)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(stream=sys.stdout)
    formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.propagate = False
    return logger



def init_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Creates a cisTopic object from fragments, consensus peaks and QC files "
            "without materialising a dense regions x cells matrix."
        )
    )
    parser.add_argument("--sample_id", type=str, metavar="<val>", help="Sample identificator")
    parser.add_argument(
        "--fragments",
        metavar="<file>",
        type=str,
        help=(
            "Specify a path to the fragments.tsv.gz file "
            "(fragments.tsv.gz should be in the same directory)"
        ),
    )
    parser.add_argument(
        "--consensus",
        metavar="<file>",
        type=str,
        help="Specify a path to the file with consensus peaks",
    )
    parser.add_argument(
        "--blacklist",
        metavar="<file>",
        type=str,
        help="Specify a path to bed file containing blacklist regions (Amemiya et al., 2019)",
    )
    parser.add_argument(
        "--qc_dir",
        metavar="<dir>",
        type=str,
        help="Specify a path to the directory with qualtiry control results",
    )
    parser.add_argument(
        "--cpus",
        metavar="<num>",
        type=int,
        help="Specify a number of cpu cores to use for interval joins",
        default=1,
    )
    parser.add_argument(
        "--unique_fragments_threshold",
        metavar="<val>",
        type=int,
        help="Threshold for number of unique fragments in peaks",
        default=None,
    )
    parser.add_argument(
        "--tss_enrichment_threshold",
        metavar="<val>",
        type=float,
        help="Threshold for TSS enrichment score",
        default=None,
    )
    parser.add_argument(
        "--frip_threshold",
        metavar="<val>",
        type=float,
        help=(
            "Threshold for fraction of reads in peaks (FRiP). "
            "If not defined the threshold will be set to 0"
        ),
        default=0,
    )
    parser.add_argument(
        "--min_frag",
        metavar="<int>",
        type=int,
        help="Minimal number of fragments in a cell for the cell to be kept. Default: 1",
        default=1,
    )
    parser.add_argument(
        "--min_cell",
        metavar="<int>",
        type=int,
        help="Minimal number of cell in which a region is detected to be kept. Default: 1",
        default=1,
    )
    parser.add_argument(
        "--is_acc",
        metavar="<int>",
        type=int,
        help="Minimal number of fragments for a region to be considered accessible. Default: 1",
        default=1,
    )
    parser.add_argument(
        "--check_for_duplicates",
        help=(
            "If no duplicate counts are provided per row in the fragments file, whether to "
            "collapse duplicates. Default: False"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--use_automatic_thresholds",
        help=(
            "Use automatic thresholds for unique fragments in peaks and TSS enrichment "
            "score as calculated by Otsu's method"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--split_pattern",
        metavar="<val>",
        type=str,
        help="Pattern to split cell barcode from sample id. Default: '___'",
        default="___",
    )
    parser.add_argument(
        "--filter_low_quality_cells",
        help="Whether to filter low quality cells based on QC metrics. Default: False",
        action="store_true",
    )
    return parser



def read_metrics(
    qc_dir: str,
    sample_id: str,
    barcodes: np.ndarray,
    filter_low_quality_cells: bool = False,
) -> DataFrame:
    fragments_stats_file = os.path.join(qc_dir, f"{sample_id}.fragments_stats_per_cb.parquet")
    if not os.path.exists(fragments_stats_file):
        raise ValueError(f'No file with path "{fragments_stats_file}" was found')

    barcodes = barcodes if barcodes.shape else [barcodes[()]]
    fragments_stats = read_parquet(fragments_stats_file, engine="pyarrow")
    fragments_stats = fragments_stats.set_index("CB")
    fragments_stats["passed_qc"] = False
    fragments_stats.loc[barcodes, "passed_qc"] = True
    if filter_low_quality_cells:
        fragments_stats = fragments_stats.loc[fragments_stats["passed_qc"]].copy()
    return fragments_stats



def _ensure_unique_strings(values: Iterable[object]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        string_value = str(value)
        if string_value in seen:
            continue
        seen.add(string_value)
        out.append(string_value)
    return out



def _region_names_to_coordinates(region_names: list[str]) -> pd.DataFrame:
    region_index = pd.Index(region_names, name="RegionIDs")
    coords = region_index.to_series(index=region_index).str.extract(
        r"^(?P<Chromosome>[^:]+):(?P<Start>\d+)-(?P<End>\d+)$"
    )
    if coords.isna().any().any():
        bad_regions = coords.index[coords.isna().any(axis=1)].tolist()
        raise ValueError(
            "The following region names do not match 'chrom:start-end': "
            + ", ".join(map(str, bad_regions[:5]))
        )
    coords["Start"] = coords["Start"].astype(np.int32)
    coords["End"] = coords["End"].astype(np.int32)
    coords.index = region_index
    return coords



def _make_binary_matrix(fragment_matrix: sparse.csr_matrix, is_acc: int) -> sparse.csr_matrix:
    binary_matrix = fragment_matrix.copy()
    if binary_matrix.nnz == 0:
        return binary_matrix
    if is_acc <= 1:
        binary_matrix.data = np.ones(binary_matrix.data.shape[0], dtype=np.int32)
    else:
        binary_matrix.data = (binary_matrix.data >= is_acc).astype(np.int32, copy=False)
        binary_matrix.eliminate_zeros()
    return binary_matrix



def _subset_list(values: list[str], idx: np.ndarray) -> list[str]:
    return [values[i] for i in idx.tolist()]



def create_cistopic_object_from_sparse_matrix(
    fragment_matrix: sparse.csr_matrix,
    cell_names: list[str],
    region_names: list[str],
    path_to_fragments: str | dict[str, str] | None,
    project: str,
    min_frag: int = 1,
    min_cell: int = 1,
    is_acc: int = 1,
    tag_cells: bool = True,
    split_pattern: str = "___",
) -> CistopicObject:
    logger = get_logger()

    if not sparse.isspmatrix_csr(fragment_matrix):
        fragment_matrix = sparse.csr_matrix(fragment_matrix, dtype=np.int32)
    else:
        fragment_matrix = fragment_matrix.astype(np.int32, copy=False)
    fragment_matrix.eliminate_zeros()

    if tag_cells:
        cell_names = [f"{cell}{split_pattern}{project}" for cell in cell_names]

    logger.info("Creating CistopicObject from sparse matrix")
    binary_matrix = _make_binary_matrix(fragment_matrix, is_acc=is_acc)

    selected_region_idx = np.flatnonzero(binary_matrix.getnnz(axis=1))
    fragment_matrix = fragment_matrix[selected_region_idx, :]
    binary_matrix = binary_matrix[selected_region_idx, :]
    region_names = _subset_list(region_names, selected_region_idx)

    selected_cell_idx = np.flatnonzero(fragment_matrix.getnnz(axis=0))
    fragment_matrix = fragment_matrix[:, selected_cell_idx]
    binary_matrix = binary_matrix[:, selected_cell_idx]
    cell_names = _subset_list(cell_names, selected_cell_idx)

    if fragment_matrix.shape[0] == 0:
        raise ValueError("No regions remain after accessibility filtering.")
    if fragment_matrix.shape[1] == 0:
        raise ValueError("No cells remain after overlap counting.")

    cisTopic_nr_frag = np.asarray(fragment_matrix.sum(axis=0)).ravel()
    cisTopic_nr_acc = np.asarray(binary_matrix.sum(axis=0)).ravel()
    cell_data = pd.DataFrame(
        {
            "cisTopic_nr_frag": cisTopic_nr_frag,
            "cisTopic_log_nr_frag": np.log10(cisTopic_nr_frag),
            "cisTopic_nr_acc": cisTopic_nr_acc,
            "cisTopic_log_nr_acc": np.log10(cisTopic_nr_acc),
            "sample_id": [project] * len(cell_names),
        },
        index=cell_names,
    )

    if min_frag != 1:
        keep_cells = cell_data["cisTopic_nr_frag"].to_numpy() >= min_frag
        fragment_matrix = fragment_matrix[:, keep_cells]
        binary_matrix = binary_matrix[:, keep_cells]
        cell_data = cell_data.loc[keep_cells].copy()
        cell_names = cell_data.index.to_list()

    if fragment_matrix.shape[1] == 0:
        raise ValueError("No cells remain after min_frag filtering.")

    region_data = _region_names_to_coordinates(region_names)
    region_data["Width"] = (region_data["End"] - region_data["Start"]).abs().astype(np.int32)
    region_data["cisTopic_nr_frag"] = np.asarray(fragment_matrix.sum(axis=1)).ravel()
    region_data["cisTopic_log_nr_frag"] = np.log10(region_data["cisTopic_nr_frag"])
    region_data["cisTopic_nr_acc"] = np.asarray(binary_matrix.sum(axis=1)).ravel()
    region_data["cisTopic_log_nr_acc"] = np.log10(region_data["cisTopic_nr_acc"])

    if min_cell != 1:
        keep_regions = region_data["cisTopic_nr_acc"].to_numpy() >= min_cell
        fragment_matrix = fragment_matrix[keep_regions, :]
        binary_matrix = binary_matrix[keep_regions, :]
        region_data = region_data.loc[keep_regions].copy()
        region_names = region_data.index.to_list()

    if fragment_matrix.shape[0] == 0:
        raise ValueError("No regions remain after min_cell filtering.")

    if path_to_fragments is None:
        path_to_fragments = {}

    cistopic_obj = CistopicObject(
        fragment_matrix=fragment_matrix,
        binary_matrix=binary_matrix,
        cell_names=cell_names,
        region_names=region_names,
        cell_data=cell_data,
        region_data=region_data,
        path_to_fragments=path_to_fragments,
        project=project,
    )
    logger.info("Done!")
    return cistopic_obj



def _load_fragments_dataframe(path_to_fragments: str, check_for_duplicates: bool) -> pd.DataFrame:
    logger = get_logger()
    logger.info("Reading fragments file")
    fragments = read_fragments_to_pyranges(
        fragments_bed_filename=path_to_fragments,
        engine="polars",
    )
    fragments_df = fragments.df

    has_score_column = "Score" in fragments_df.columns
    fragments_df = fragments_df.loc[:, ["Chromosome", "Start", "End", "Name"]].copy()

    if (not has_score_column) and check_for_duplicates:
        logger.info("Collapsing duplicate fragment rows")
        fragments_df = fragments_df.drop_duplicates(
            subset=["Chromosome", "Start", "End", "Name"],
            keep="first",
            ignore_index=True,
        )

    return fragments_df



def _prepare_candidate_barcodes(
    fragments_df: pd.DataFrame,
    metrics: str | pd.DataFrame | None,
    valid_bc: list[str] | None,
) -> tuple[pd.DataFrame | None, list[str]]:
    logger = get_logger()

    if metrics is not None:
        logger.info("metrics provided!")
        if isinstance(metrics, str):
            metrics = pd.read_csv(metrics)
        else:
            metrics = metrics.copy()

        if "is__cell_barcode" in metrics.columns:
            metrics = metrics.loc[metrics["is__cell_barcode"] == 1].copy()
            metrics.index = metrics["barcode"].astype(str)
            metrics = metrics.iloc[:, 2:]
        else:
            metrics.index = metrics.index.astype(str)

        metrics = metrics.loc[~metrics.index.duplicated(keep="first")].copy()
        candidate_barcodes = _ensure_unique_strings(metrics.index)
        return metrics, candidate_barcodes

    if valid_bc is not None:
        candidate_barcodes = _ensure_unique_strings(valid_bc)
        return None, candidate_barcodes

    candidate_barcodes = _ensure_unique_strings(fragments_df["Name"].tolist())
    return None, candidate_barcodes



def _filter_fragments_to_candidate_barcodes(
    fragments_df: pd.DataFrame,
    candidate_barcodes: list[str],
) -> pd.DataFrame:
    logger = get_logger()
    logger.info("Filtering fragments to candidate barcodes")
    fragments_df = fragments_df.copy()
    fragments_df["Name"] = pd.Categorical(
        values=fragments_df["Name"],
        categories=candidate_barcodes,
        ordered=True,
    )
    fragments_df["CellIdx"] = fragments_df["Name"].cat.codes.astype(np.int32)
    fragments_df = fragments_df.loc[fragments_df["CellIdx"] >= 0, ["Chromosome", "Start", "End", "Name", "CellIdx"]].copy()
    return fragments_df



def _make_unique_fragment_counts_per_barcode(fragments_df: pd.DataFrame) -> pd.DataFrame:
    counts = fragments_df.groupby("Name", sort=False, observed=True).size().astype(np.int32)
    counts.index = counts.index.astype(str)
    counts_df = pd.DataFrame({"Unique_nr_frag": counts})
    counts_df["barcode"] = counts_df.index.to_list()
    return counts_df



def _load_regions_dataframe(
    path_to_regions: str,
    path_to_blacklist: str | None,
) -> pd.DataFrame:
    logger = get_logger()
    logger.info("Reading consensus regions")
    regions = pr.read_bed(path_to_regions)
    regions_df = regions.df.loc[:, ["Chromosome", "Start", "End"]].copy()
    regions_df["regionID"] = (
        regions_df["Chromosome"].astype(str)
        + ":"
        + regions_df["Start"].astype(str)
        + "-"
        + regions_df["End"].astype(str)
    )

    if isinstance(path_to_blacklist, str):
        logger.info("Removing blacklisted regions before counting")
        blacklist = pr.read_bed(path_to_blacklist)
        regions = pr.PyRanges(regions_df).overlap(blacklist, invert=True)
        regions_df = regions.df.loc[:, ["Chromosome", "Start", "End", "regionID"]].copy()

    return regions_df



def _build_sparse_fragment_matrix_by_chromosome(
    fragments_df: pd.DataFrame,
    regions_df: pd.DataFrame,
    n_cells: int,
    n_cpu: int,
) -> tuple[sparse.csr_matrix, list[str]]:
    logger = get_logger()
    blocks: list[sparse.csr_matrix] = []
    kept_region_names: list[str] = []

    region_chromosomes = pd.Index(regions_df["Chromosome"]).drop_duplicates().tolist()

    for chromosome in region_chromosomes:
        regions_chr = regions_df.loc[regions_df["Chromosome"] == chromosome].copy()
        fragments_chr = fragments_df.loc[fragments_df["Chromosome"] == chromosome, ["Chromosome", "Start", "End", "CellIdx"]]

        if regions_chr.empty or fragments_chr.empty:
            continue

        regions_chr["LocalRow"] = np.arange(regions_chr.shape[0], dtype=np.int32)

        joined = pr.PyRanges(regions_chr.loc[:, ["Chromosome", "Start", "End", "LocalRow"]]).join(
            pr.PyRanges(fragments_chr),
            nb_cpu=max(1, n_cpu),
        )
        joined_df = joined.df

        if joined_df.empty:
            logger.info("Chromosome %s: no overlaps", chromosome)
            continue

        row = joined_df["LocalRow"].to_numpy(dtype=np.int64, copy=False)
        col = joined_df["CellIdx"].to_numpy(dtype=np.int64, copy=False)
        data = np.ones(row.shape[0], dtype=np.int32)

        block = sparse.coo_matrix(
            (data, (row, col)),
            shape=(regions_chr.shape[0], n_cells),
            dtype=np.int32,
        ).tocsr()
        block.eliminate_zeros()

        nonzero_rows = np.flatnonzero(block.getnnz(axis=1))
        if nonzero_rows.size == 0:
            logger.info("Chromosome %s: overlaps collapsed to zero signal", chromosome)
            continue

        block = block[nonzero_rows, :]
        blocks.append(block)
        kept_region_names.extend(regions_chr.iloc[nonzero_rows]["regionID"].astype(str).tolist())

        logger.info(
            "Chromosome %s: %d overlaps -> %d non-empty regions",
            chromosome,
            row.shape[0],
            nonzero_rows.size,
        )

        del joined, joined_df, row, col, data, regions_chr, fragments_chr, nonzero_rows
        gc.collect()

    if not blocks:
        raise ValueError("No fragment overlaps were found for the provided regions and barcodes.")

    fragment_matrix = sparse.vstack(blocks, format="csr", dtype=np.int32)
    fragment_matrix.sort_indices()
    logger.info(
        "Sparse fragment matrix shape: %d regions x %d cells with %d non-zero entries",
        fragment_matrix.shape[0],
        fragment_matrix.shape[1],
        fragment_matrix.nnz,
    )
    return fragment_matrix, kept_region_names



def create_cistopic_object_from_fragments_memory_optimized(
    path_to_fragments: str,
    path_to_regions: str,
    path_to_blacklist: str | None = None,
    metrics: str | pd.DataFrame | None = None,
    valid_bc: list[str] | None = None,
    n_cpu: int = 1,
    min_frag: int = 1,
    min_cell: int = 1,
    is_acc: int = 1,
    check_for_duplicates: bool = True,
    project: str = "cisTopic",
    split_pattern: str = "___",
) -> tuple[CistopicObject, pd.DataFrame | None]:
    logger = get_logger()
    logger.info("Reading data for %s", project)

    fragments_df = _load_fragments_dataframe(
        path_to_fragments=path_to_fragments,
        check_for_duplicates=check_for_duplicates,
    )

    metrics_df, candidate_barcodes = _prepare_candidate_barcodes(
        fragments_df=fragments_df,
        metrics=metrics,
        valid_bc=valid_bc,
    )
    if not candidate_barcodes:
        raise ValueError("No candidate barcodes were available after QC filtering.")
    logger.info("Candidate barcodes retained: %d", len(candidate_barcodes))

    fragments_df = _filter_fragments_to_candidate_barcodes(
        fragments_df=fragments_df,
        candidate_barcodes=candidate_barcodes,
    )
    if fragments_df.empty:
        raise ValueError("No fragment rows remain after barcode filtering.")
    logger.info("Fragment rows retained after barcode filtering: %d", fragments_df.shape[0])

    unique_fragment_counts = None
    if metrics_df is None:
        logger.info("Counting number of unique fragments per barcode (Unique_nr_frag)")
        unique_fragment_counts = _make_unique_fragment_counts_per_barcode(fragments_df)

    regions_df = _load_regions_dataframe(
        path_to_regions=path_to_regions,
        path_to_blacklist=path_to_blacklist,
    )

    logger.info("Building sparse fragment matrix chromosome-by-chromosome")
    fragment_matrix, region_names = _build_sparse_fragment_matrix_by_chromosome(
        fragments_df=fragments_df.loc[:, ["Chromosome", "Start", "End", "CellIdx"]],
        regions_df=regions_df,
        n_cells=len(candidate_barcodes),
        n_cpu=n_cpu,
    )
    del fragments_df, regions_df
    gc.collect()

    cistopic_obj = create_cistopic_object_from_sparse_matrix(
        fragment_matrix=fragment_matrix,
        cell_names=candidate_barcodes,
        region_names=region_names,
        path_to_fragments={project: path_to_fragments},
        project=project,
        min_frag=min_frag,
        min_cell=min_cell,
        is_acc=is_acc,
        tag_cells=True,
        split_pattern=split_pattern,
    )

    if metrics_df is not None:
        metrics_df = metrics_df.copy()
        metrics_df["barcode"] = metrics_df.index.to_list()
        cistopic_obj.add_cell_data(metrics_df, split_pattern)
    elif unique_fragment_counts is not None:
        cistopic_obj.add_cell_data(unique_fragment_counts, split_pattern)

    return cistopic_obj, metrics_df if metrics_df is not None else unique_fragment_counts



def main() -> None:
    parser = init_parser()
    args = parser.parse_args()

    barcodes, thresholds = get_barcodes_passing_qc_for_sample(
        sample_id=args.sample_id,
        pycistopic_qc_output_dir=args.qc_dir,
        unique_fragments_threshold=args.unique_fragments_threshold,
        tss_enrichment_threshold=args.tss_enrichment_threshold,
        frip_threshold=args.frip_threshold,
        use_automatic_thresholds=args.use_automatic_thresholds,
    )

    _, otsu_thresholds = get_barcodes_passing_qc_for_sample(
        sample_id=args.sample_id,
        pycistopic_qc_output_dir=args.qc_dir,
        use_automatic_thresholds=args.use_automatic_thresholds,
    )

    fragments_stats = read_metrics(
        args.qc_dir,
        args.sample_id,
        barcodes,
        filter_low_quality_cells=args.filter_low_quality_cells,
    )

    cistopic_obj, _ = create_cistopic_object_from_fragments_memory_optimized(
        path_to_fragments=args.fragments,
        path_to_regions=args.consensus,
        path_to_blacklist=args.blacklist,
        metrics=fragments_stats,
        valid_bc=None,
        n_cpu=args.cpus or 1,
        min_frag=args.min_frag,
        min_cell=args.min_cell,
        is_acc=args.is_acc,
        check_for_duplicates=args.check_for_duplicates,
        project=args.sample_id,
        split_pattern=args.split_pattern,
    )

    with open(f"{args.sample_id}_cistopic_obj.pkl", "wb") as handle:
        pickle.dump(cistopic_obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(f"{args.sample_id}_good_cells.txt", "w") as handle:
        handle.write("\n".join(barcodes.tolist()))

    with open(f"{args.sample_id}_thresholds.json", "w") as handle:
        json.dump(thresholds, handle)

    adata = ad.AnnData(
        X=cistopic_obj.fragment_matrix.T,
        obs=cistopic_obj.cell_data.infer_objects(),
        var=cistopic_obj.region_data.infer_objects(),
        layers={"binary": cistopic_obj.binary_matrix.T},
        uns={"qc_thresholds": thresholds, "otsu_thresholds": otsu_thresholds},
    )
    adata.write_h5ad(f"{args.sample_id}.h5ad")


if __name__ == "__main__":
    main()