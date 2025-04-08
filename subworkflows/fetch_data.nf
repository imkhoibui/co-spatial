include { PREFETCH                      } from "${projectDir}/modules/local/prefetch/main.nf"
include { FASTERQ_DUMP                  } from "${projectDir}/modules/local/fasterq-dump/main.nf"

workflow FETCH_DATA {
    take:
        ch_input
        ch_fastq_out

    main:

        ch_input.combine(ch_fastq_out)
            .map { meta, fastq_out -> 
                def meta_id = meta["sample_id"]
                def fastq_path = [fastq_out.toString(), meta["experiment"].toString(), meta_id].join("/")
                return [meta: meta_id, fastq_out: fastq_path]
            }
            .set { ch_input_fasterq_dump }

        PREFETCH(
            ch_input
        )

        PREFETCH.out.sra.view()


        // FASTERQ_DUMP(
        //     ch_fastq_out
        // )
}