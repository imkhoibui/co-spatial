include { FETCH_DATA                             } from "${projectDir}/subworkflows/local/fetch_data.nf"
include { SPATIAL_RNA                            } from "${projectDir}/subworkflows/local/spatial_rna.nf"
include { SPATIAL_ATAC                           } from "${projectDir}/subworkflows/local/spatial_atac.nf"

workflow SPATIAL_ANALYSIS {
    ch_input                = Channel.fromPath(params.input).splitCsv( header: true )
    ch_fastq_out            = Channel.fromPath(params.fastq_out)
    ch_genome_meta          = Channel.of(params.meta)

    ch_spatial_barcodes     = Channel.fromPath(params.spatial_barcodes)
    ch_ref_map              = ch_genome_meta.combine(Channel.fromPath(params.ref_map))
    ch_ref_annotation       = Channel.fromPath(params.ref_annotation)
    ch_ref_atac_genome      = Channel.fromPath(params.ref_atac_genome)
    ch_genome_fasta         = Channel.fromPath(params.whole_genome_fasta)
    ch_genome_fasta_files   = ch_genome_meta.combine(ch_genome_fasta)

    // Module to fetch SRA data
    if ( !params.skip_fetch_data ) {
        FETCH_DATA(
            ch_input,
            ch_fastq_out
        )
    }

    ch_input.combine(ch_fastq_out)
        .map { meta, fastq_out -> 
            def meta_id = meta["sample_id"]
            def experiment = meta["experiment"]
            def fastq_path = [fastq_out.toString(), meta_id].join("/")
            return [meta_id, experiment, fastq_path]
        }
        .branch {
            atacseq: it[1] == "atacseq"
            rnaseq: it[1] == "rnaseq"
        }
        .set { ch_input_fastq }

    ch_input_atacseq = ch_input_fastq.atacseq
        .map { meta_id, experiment, fastq_path ->
            def fastq1_path = file([fastq_path, meta_id + "_1.fastq.gz"].join("/"), checkIfExists: true)
            def fastq2_path = file([fastq_path, meta_id + "_2.fastq.gz"].join("/"), checkIfExists: true)
            return [meta_id, fastq1_path, fastq2_path]
        }
    ch_input_rnaseq = ch_input_fastq.rnaseq
        .map { meta_id, experiment, fastq_path ->
            def fastq1_path = file([fastq_path, meta_id + "_1.fastq.gz"].join("/"), checkIfExists: true)
            def fastq2_path = file([fastq_path, meta_id + "_2.fastq.gz"].join("/"), checkIfExists: true)
            return [meta_id, fastq1_path, fastq2_path]
        }

    // Running spatial RNAseq module
    SPATIAL_RNA(
        ch_input_rnaseq,
        ch_spatial_barcodes,
        ch_ref_map,
        ch_ref_annotation,
        ch_genome_fasta_files
    )

    // Running spatial ATACseq module
    SPATIAL_ATAC(
        ch_input_atacseq,
        ch_ref_atac_genome
    )
}