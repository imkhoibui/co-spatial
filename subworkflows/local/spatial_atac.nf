include { ATAC_PREPROCESS_FILTER_PRIMER                 } from "${projectDir}/modules/local/filter_primer.nf"
include { ATAC_PREPROCESS_RENAME                        } from "${projectDir}/modules/local/atac_preprocess.nf"

workflow SPATIAL_ATAC {
    take:
        ch_input_actacseq
        ch_genome_ref

    main:

        ATAC_PREPROCESS_FILTER_PRIMER(
            ch_input_actacseq
        )

        ATAC_PREPROCESS_RENAME(
            ATAC_PREPROCESS_FILTER_PRIMER.out.filter
        )



}