include { SPATIAL_PROCESS } from "${projectDir}/modules/local/spatial_process.nf"

workflow SPATIAL_RNA_ATAC {
    take:
        ch_st_pipeline
        ch_atac_outputs
        ch_tissue_dir 
        ch_spatial_barcodes

    main:
        SPATIAL_PROCESS(
            ch_st_pipeline,
            ch_atac_outputs,
            ch_tissue_dir,
            ch_spatial_barcodes
        )

    emit:
        SPATIAL_PROCESS.out.rds
}