include { FETCH_DATA                             } from "${projectDir}/subworkflows/fetch_data.nf"

workflow SPATIAL_ANALYSIS {
    ch_input            = Channel.fromPath(params.input).splitCsv( header: true )
    ch_fastq_out        = Channel.fromPath(params.fastq_out)
    ch_outdir           = Channel.fromPath(params.outdir)

    FETCH_DATA(
        ch_input,
        ch_fastq_out
    )

}