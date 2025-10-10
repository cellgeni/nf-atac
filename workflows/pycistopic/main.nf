include { CISTOPIC_PEAKCALLING } from '../../subworkflows/local/cistopic_peakcalling'
include { CISTOPIC_INFERPEAKS } from '../../subworkflows/local/cistopic_inferpeaks'
include { CISTOPIC_ATTACHGEX } from '../../subworkflows/local/cistopic_attachgex'


workflow  PYCISTOPIC {
    take:
    sample_table
    celltypes
    pseudobulk_peaks
    atac_adata
    chromsizes
    blacklist
    tss_bed
    callPeaksFlag
    inferConsensusFlag
    attachGEXFlag
    
    main:
    versions         = Channel.empty()
    pseudobulk_peaks = pseudobulk_peaks.splitCsv(header: true, sep:',').map { row -> tuple( [id: row.celltype, fragments: row.fragments], file( row.path, checkIfExists: true ) ) }
    atac_adata       = atac_adata.splitCsv(header: true, sep:'\t').map { row -> tuple( [id: row.sample_id], file( row.path ) ) }

    // Create pseudobulk, call peaks and update sample table
    if ( callPeaksFlag ) {
        CISTOPIC_PEAKCALLING(
            sample_table,
            celltypes,
            chromsizes
        )

        versions = versions.mix(CISTOPIC_PEAKCALLING.out.versions)
        // get pseudobulk peaks and updated sample table
        peak_metadata = CISTOPIC_PEAKCALLING.out.peaks
        updated_samples = CISTOPIC_PEAKCALLING.out.sample_table
    }

    if ( inferConsensusFlag ) {

        if ( ! callPeaksFlag ) {
            updated_samples = sample_table
            peak_metadata = pseudobulk_peaks
        }

        CISTOPIC_INFERPEAKS(
            peak_metadata,
            updated_samples,
            chromsizes,
            blacklist,
            tss_bed
        )
        versions = versions.mix(CISTOPIC_INFERPEAKS.out.versions)
    }

    if ( attachGEXFlag ) {
        if ( ! callPeaksFlag ) {
            updated_samples = sample_table
        }

        CISTOPIC_ATTACHGEX(
            updated_samples,
            celltypes,
            CISTOPIC_INFERPEAKS.out.anndata
        )

        versions = versions.mix(CISTOPIC_ATTACHGEX.out.versions)
    }

    emit:
    pseudobulk   = callPeaksFlag ? CISTOPIC_PEAKCALLING.out.pseudobulk : Channel.empty()
    peaks        = callPeaksFlag ? CISTOPIC_PEAKCALLING.out.peaks : Channel.empty()
    consensus    = inferConsensusFlag ? CISTOPIC_INFERPEAKS.out.consensus : Channel.empty()
    cistopic     = inferConsensusFlag ? CISTOPIC_INFERPEAKS.out.cistopic : Channel.empty()
    atac_anndata = inferConsensusFlag ? CISTOPIC_INFERPEAKS.out.anndata : Channel.empty()
    coupled_h5mu = attachGEXFlag ? CISTOPIC_ATTACHGEX.out.h5mu : Channel.empty()
    coupled_h5ad = attachGEXFlag ? CISTOPIC_ATTACHGEX.out.h5ad : Channel.empty()
    versions     = versions
}