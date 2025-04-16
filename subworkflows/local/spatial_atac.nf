include { CUTADAPT as CUTADAPT_PRIMER                   } from "${projectDir}/modules/local/cutadapt/main.nf"
include { CUTADAPT as CUTADAPT_LINKER1                  } from "${projectDir}/modules/local/cutadapt/main.nf"
include { CUTADAPT as CUTADAPT_LINKER2                  } from "${projectDir}/modules/local/cutadapt/main.nf"

include { ATAC_PREPROCESS                               } from "${projectDir}/modules/local/atac_preprocess.nf"
include { CELLRANGER_MKREF                              } from "${projectDir}/modules/local/cellranger/mkref.nf"
include { CELLRANGER                                    } from "${projectDir}/modules/local/cellranger/main.nf"

workflow SPATIAL_ATAC {
    take:
        ch_input_actacseq
        ch_ref_atac_genome

    main:

        // Trim primer & linker sequences
        CUTADAPT_PRIMER(
            ch_input_actacseq
        )
        CUTADAPT_LINKER1(
            CUTADAPT_PRIMER.out.fastq
        )
        CUTADAPT_LINKER2(
            CUTADAPT_LINKER1.out.fastq
        )

        // Perform preprocess & adapt for cellranger ATAC format
        ATAC_PREPROCESS(
            CUTADAPT_LINKER2.out.fastq
        )
        ch_cellranger_input = ATAC_PREPROCESS.out.fastq
        
        // Running CellRanger ATAC
        if (!params.skip_make_ref) {
            CELLRANGER_MKREF()
            ch_ref_atac_genome = CELLRANGER_MKREF().out.build
        }

        CELLRANGER(
            ch_cellranger_input,
            ch_ref_atac_genome
        )

    emit: 
        atac_outputs = CELLRANGER.out.outputs
}