include { CISTOPIC_PEAKCALLING } from '../../subworkflows/local/cistopic_peakcalling'
include { CISTOPIC_INFERPEAKS } from '../../subworkflows/local/cistopic_inferpeaks'


workflow  PYCISTOPIC {
    take:
    sample_table
    celltypes
    chromsizes
    blacklist
    tss_bed
    callPeaksFlag
    inferConsensusFlag
    fragments_filename
    narrowPeaks_dir
    
    main:
    versions = Channel.empty()

    // Create pseudobulk, call peaks and update sample table
    if ( callPeaksFlag ) {
        CISTOPIC_PEAKCALLING(
            sample_table,
            celltypes,
            chromsizes,
            fragments_filename,
            narrowPeaks_dir
        )

        versions = versions.mix(CISTOPIC_PEAKCALLING.out.versions)
        // get pseudobulk peaks and updated sample table
        peak_metadata = CISTOPIC_PEAKCALLING.out.peaks
        sample_table = CISTOPIC_PEAKCALLING.out.sample_table
    } else {
        sample_table = Channel.fromPath(sample_table, checkIfExists: true)
        peak_metadata = Channel.fromPath(celltypes, checkIfExists: true).splitCsv(skip:1, sep:'\t')
    }

    if ( inferConsensusFlag ) {
        CISTOPIC_INFERPEAKS(
            peak_metadata,
            sample_table,
            chromsizes,
            blacklist,
            tss_bed,
            fragments_filename
        )
        versions = versions.mix(CISTOPIC_INFERPEAKS.out.versions)
    }

    emit:
    pseudobulk = callPeaksFlag ? CISTOPIC_PEAKCALLING.out.pseudobulk : Channel.empty()
    peaks      = callPeaksFlag ? CISTOPIC_PEAKCALLING.out.peaks : Channel.empty()
    consensus  = inferConsensusFlag ? CISTOPIC_INFERPEAKS.out.consensus : Channel.empty()
    cistopic   = inferConsensusFlag ? CISTOPIC_INFERPEAKS.out.cistopic : Channel.empty()
    anndata    = inferConsensusFlag ? CISTOPIC_INFERPEAKS.out.anndata : Channel.empty()
    versions   = versions
}