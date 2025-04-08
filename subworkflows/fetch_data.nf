include { PREFETCH                      } from "${projectDir}/modules/local/prefetch/main.nf"
include { FASTERQ_DUMP                  } from "${projectDir}/modules/local/fasterq-dump/main.nf"

workflow FETCH_DATA {
    take:
        ch_input
        ch_fastq_out

    main:
        ch_asc_id = ch_input.map { 
            [meta: it[0].take(3), asc: it[0]] 
        }
        ch_fastq_out = ch_fastq_out + "/" + ch_input.map { it[1] }

        ch_asc_id.view()
        ch_fastq_out.view()

        PREFETCH(
            ch_asc_id
        )

        FASTERQ_DUMP(
            ch_fastq_out
        )
}