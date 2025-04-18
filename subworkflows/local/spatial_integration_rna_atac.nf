include { SIGNAC_RNA } from "${projectDir}/modules/local/signac/signac_rna.nf"
include { ARCH_R } from "${projectDir}/modules/local/archR/main.nf"

workflow SPATIAL_RNA_ATAC {
    take:
        ch_st_pipeline
        ch_atac_outputs
        ch_tissue_dir 
        ch_spatial_barcodes

    main:
        ch_st_pipeline = ch_st_pipeline
        SIGNAC_RNA(
            ch_st_pipeline,
            ch_st_tissue_dir,
            ch_spatial_barcodes
        )

        // ARCH_R(
        //     ch_st_pipeline,
        //     ch_atac_outputs
        // )

    emit:
        outputs = ch_st_pipeline
}