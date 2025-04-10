include { FETCH_DATA                             } from "${projectDir}/subworkflows/local/fetch_data.nf"
include { SPATIAL_RNA                            } from "${projectDir}/subworkflows/local/spatial_rna.nf"

workflow SPATIAL_ANALYSIS {
    ch_input            = Channel.fromPath(params.input).splitCsv( header: true )
    ch_fastq_out        = Channel.fromPath(params.fastq_out)
    ch_outdir           = Channel.fromPath(params.outdir)

    if ( !params.skip_fetch_data ) {
        FETCH_DATA(
            ch_input,
            ch_fastq_out
        )
    }

    ch_input.combine(ch_fastq_out)
        .map { meta, fastq_out -> 
            def meta_id = meta["sample_id"]
            def fastq_path = [fastq_out.toString(), meta["experiment"].toString(), meta_id].join("/")
            return [meta_id, fastq_path]
        }
        .set { ch_input_fastq }

    ch_input_fastq.view()

    // SCRNA_PREPROCESS(
    //     ch_input_fastq
    // )


}