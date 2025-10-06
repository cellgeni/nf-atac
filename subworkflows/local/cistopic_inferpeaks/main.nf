include { CISTOPIC_INFERCONSENSUS } from '../../../modules/local/cistopic/inferconsensus'
include { CISTOPIC_QUALITYCONTROL } from '../../../modules/local/cistopic/qualitycontrol'
include { CISTOPIC_CREATEOBJECT } from '../../../modules/local/cistopic/createobject'
include { CISTOPIC_COMBINEOBJECTS } from '../../../modules/local/cistopic/combineobjects'
include { ANNDATA_CONCAT } from '../../../modules/local/anndata/concat'


workflow CISTOPIC_INFERPEAKS {
    take:
        peaks
        sample_table
        chromsizes
        blacklist
        tss_bed
    main:
        // STEP 0: Prepare inputs
        // Get fragment paths from sample table
        fragments = sample_table
            .splitCsv(header: true)
            .map{ meta -> 
                [
                    [id: meta.sample_id, fragments: meta.fragments],
                    file( "${meta.filedir}/*fragments.tsv.gz" )[0],
                    file( "${meta.filedir}/*fragments.tsv.gz.tbi" )[0]
                ]
            }

        // Get peak paths from peaks channel
        narrowPeaks = peaks.toSortedList()
            .map{ list -> 
                def celltype_names = list.collect{ it[0].id }
                def narrowpeak_files = list.collect{ it[1] }
                return [ celltype_names, narrowpeak_files ]
            }
        
        // STEP 1: Get consensus peaks
        CISTOPIC_INFERCONSENSUS(
            narrowPeaks,
            chromsizes,
            blacklist
        )

        consensus = CISTOPIC_INFERCONSENSUS.out.bed.collect()

        // STEP 2: Perform QC
        CISTOPIC_QUALITYCONTROL(fragments, consensus, tss_bed)

        // SPET 3: Create cisTopic object
        CISTOPIC_CREATEOBJECT(
            CISTOPIC_QUALITYCONTROL.out.qc,
            consensus,
            blacklist
        )
        
        cistopic_objects = CISTOPIC_CREATEOBJECT.out.pkl
            .toSortedList()
            .branch {
                it ->
                combine_objects: it.size() > 1
                sample: true
            }
        
        anndata_objects = CISTOPIC_CREATEOBJECT.out.h5ad
            .toSortedList()
            .branch {
                it ->
                combine_objects: it.size() > 1
                sample: true
            }

        // STEP 4: Combine cisTopic objects if needed
        CISTOPIC_COMBINEOBJECTS(cistopic_objects.combine_objects.transpose().toList())
        ANNDATA_CONCAT(anndata_objects.combine_objects.transpose().toList())

        // STEP 5: Collect versions
        versions = CISTOPIC_INFERCONSENSUS.out.versions
            .mix(
                CISTOPIC_QUALITYCONTROL.out.versions.first(),
                CISTOPIC_CREATEOBJECT.out.versions.first(),
                CISTOPIC_COMBINEOBJECTS.out.versions.ifEmpty{ Channel.empty() },
                ANNDATA_CONCAT.out.versions.ifEmpty{ Channel.empty() }
            )

        emit:
        consensus = consensus
        cistopic  = cistopic_objects.sample.mix(CISTOPIC_COMBINEOBJECTS.out.pkl)
        anndata   = anndata_objects.sample.mix(ANNDATA_CONCAT.out.h5ad)
        versions  = versions

}