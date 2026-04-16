include { CISTOPIC_INFERCONSENSUS } from '../../../modules/local/cistopic/inferconsensus'
include { CISTOPIC_QUALITYCONTROL } from '../../../modules/local/cistopic/qualitycontrol'
include { CISTOPIC_CREATEOBJECT } from '../../../modules/local/cistopic/createobject'
include { CISTOPIC_COMBINEOBJECTS } from '../../../modules/local/cistopic/combineobjects'
include { ANNDATA_ATTACHCELLTYPES as ANNDATA_ATTACHCELLTYPES_ATAC } from '../../../modules/local/anndata/attachcelltypes'
include { ANNDATA_CONCAT } from '../../../modules/local/anndata/concat'


workflow CISTOPIC_INFERPEAKS {
    take:
        peaks
        chromsizes
        blacklist
    main:

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

        // STEP 2: Collect versions
        versions = CISTOPIC_INFERCONSENSUS.out.versions
        emit:
        consensus       = consensus
        versions        = versions

}