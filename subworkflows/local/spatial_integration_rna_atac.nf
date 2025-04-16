workflow SPATIAL_RNA_ATAC {
    take:
        ch_st_pipeline
        ch_atac_outputs

    main:
        ch_st_pipeline = ch_atac_outputs

    emit:
        outputs = ch_st_pipeline
}