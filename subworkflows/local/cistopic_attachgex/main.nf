include { ANNDATA_TOH5AD } from '../../../modules/local/anndata/toh5ad'
include { ANNDATA_ATTACHCELLTYPES as ANNDATA_ATTACHCELLTYPES_GEX } from '../../../modules/local/anndata/attachcelltypes'
include { ANNDATA_COUPLEMULTIOME } from '../../../modules/local/anndata/couplemultiome'

workflow CISTOPIC_ATTACHGEX {
    take:
    sample_table
    celltypes
    atac_adata

    main:
    // STEP 0: Get mtx paths from sample table
    mtx_dirs = sample_table
        .splitCsv(header: true)
        .map{ row ->
            def mtx_path = file( "${row.path}/raw_feature_bc_matrix" )
            def meta = [id: row.sample_id, fragments: row.fragments ? row.fragments : null]
            if ( ! mtx_path.exists() ) {
                error "Matrix directory not found: ${mtx_path}. Please check your sample_table."
            } 
            tuple( meta, mtx_path )
        }
    
    // STEP 1: Create anndata for GEX
    ANNDATA_TOH5AD(mtx_dirs)

    // STEP 2: Attach celltype annotation
    ANNDATA_ATTACHCELLTYPES_GEX(ANNDATA_TOH5AD.out.h5ad, celltypes)

    // STEP 3: Attach GEX to ATAC

    // Prepare the channels
    atac_adata = atac_adata.map { meta, path -> tuple( [id: meta.id], path ) }
    gex_adata = ANNDATA_ATTACHCELLTYPES_GEX.out.h5ad.map { meta, path -> tuple( [id: meta.id], path ) }
    gex_atac_pairs = atac_adata.join(gex_adata, failOnMismatch: true)
    
    // Attach GEX to ATAC
    ANNDATA_COUPLEMULTIOME(gex_atac_pairs)

    // STEP 4: Collect files
    versions = ANNDATA_TOH5AD.out.versions.first().mix(
        ANNDATA_ATTACHCELLTYPES_GEX.out.versions.first(),
        ANNDATA_COUPLEMULTIOME.out.versions.first()
        )

    emit:
    versions = versions
    h5mu = ANNDATA_COUPLEMULTIOME.out.h5mu
    h5ad = ANNDATA_COUPLEMULTIOME.out.h5ad

    



}