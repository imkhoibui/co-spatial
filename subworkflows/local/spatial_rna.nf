include { STAR_GENOMEGENERATE       } from "${projectDir}/modules/local/star/genome_generate.nf"

include { DBIT_PREPROCESS           } from "${projectDir}/modules/local/dbit_preprocess.nf"
include { ST_PIPELINE               } from "${projectDir}/modules/local/st-pipeline/main.nf"

workflow SPATIAL_RNA {
    take:
        ch_input_rnaseq
        ch_spatial_barcodes
        ch_ref_map
        ch_ref_annotation
        ch_genome_fasta_files

    main:
    DBIT_PREPROCESS{
        ch_input_rnaseq
    }

    if ( !params.skip_star_genome ) {
        STAR_GENOMEGENERATE(
            ch_genome_fasta_files,
            ch_ref_annotation
        )

        ch_ref_map = STAR_GENOMEGENERATE.out.index
    }

    ST_PIPELINE(
        DBIT_PREPROCESS.out.fastq_processed,
        ch_spatial_barcodes,
        ch_ref_map,
        ch_ref_annotation
    )

    emit:
        st_pipeline = ST_PIPELINE.out.st_pipeline
}