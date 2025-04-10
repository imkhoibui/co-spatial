include { ST_PIPELINE               } from "${projectDir}/modules/local/st-pipeline/main.nf"

workflow SPATIAL_RNA {
    take:
        ch_fastq_rna
        ch_genome_index
        ch_gff


    ST_PIPELINE(
        ch_fastq_rna,
        ch_genome_index,
        ch_gff,
        ch_exp_name
    )
}