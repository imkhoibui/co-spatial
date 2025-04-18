# CO-SPATIAL: A Nextflow pipeline to perform spatial-ATAC-RNA seq data integration

Co-spatial is a Nextflow pipeline, inspired by [nf-core](https://nf-co.re) practices, to generate spatial
ATAC-seq and spatial RNA-seq-ready counts & peaks outputs for joint spatial ATAC-RNA sequencing downstream analysis workflows.

With Nextflow, user can simultaneously process both spatial-ATACseq and spatial-RNAseq data from the same sample, this method is made possible through (Deterministic co-Barcoding in Tissue), where single cells can now contain information from both spatial transcriptomics and spatial-ATAC-seq or spatial-CUT&Tag-seq.

### Requirements
[Nextflow](https://www.nextflow.io/docs/stable/install.html)

[Docker](https://docs.docker.com/engine/install/)

This pipeline requires high computing resources to run efficiently, please modify `base.config` if needed.

### Workflows:
1. By providing a SRA accession, the data is `prefetch` and `fasterq-dump` into FASTQ files from spatial-omics experiments.
2. The downloaded files are "branched" into either spatial-RNA or spatial-ATAC for its separate processing.
    2a. For sp-RNAseq data, FASTQ files are preprocessed to extract barcode A, barcode B & UMI (read2). They are then processed using [ST_pipeline](https://github.com/jfnavarro/st_pipeline), where it is aligned to [STAR](https://github.com/alexdobin/STAR) using a reference genome & annotation.
    2b. For sp-ATACseq data, FASTQ files are preprocessed and aligned to a reference genome. The count data is then converted into 10X Genomics CellRanger-ATAC format which then can be counted using [cellranger-atac](https://github.com/10XGenomics/cellranger-atac)

### Inputs
In order to successfully run the pipeline, you must provide a `samplesheet.csv` file with the following content:

```
id,name,experiment
SRR22561636,ME13_50um,RNA
SRR22565186,ME13_50um,ATAC
```

### Running the pipeline
The usual run would have the following format:

```
nextflow run main.nf \
    -profile docker \
    --input <path/to/samplesheet.csv> \
    --fastq_out <path/to/fastq-output-dir> \
    --spatial_barcodes <path/to/spatial_barcodes.txt> \
    --ref_map <path/to/reference-genome-dir> \
    --ref_annotation <path/to/references/genes.gtf> \
    --ref_atac_genome <path/to/references/cellranger-atac> \
    --outdir <output-dir> 
```

For convenience of testing, I have created a test profile to run the pipeline on a very small set of data (10-mil reads) of 2 samples:

```
nextflow run main.nf \
    -profile test_subdata
```

### Future updates
- Integrate a CUT&Tag-RNA sequencing subworkflow
- Incorporate AWSBatch / HPC option (if I have access to one)

### Citations
This pipeline is inspired by the publication: *Spatial epigenome-transcriptome co-profiling of mammalian tissues*

> **Spatial epigenome-transcriptome co-profiling of mammalian tissues.**
>
> Di Zhang, Yanxiang Deng, Petra Kukanja, Eneritz Agirre, Marek Bartosovic, Mingze Dong, Cong Ma, Sai Ma, Graham Su, Shuozhen Bao, Yang Liu, Yang Xiao, Gorazd B. Rosoklija, Andrew J. Dwork, J. John Mann, Kam W. Leong, Maura Boldrini, Liya Wang, Maximilian Haeussler, Benjamin J. Raphael, Yuval Kluger, Gonçalo Castelo-Branco & Rong Fan 
>
> _Nature_ 2023 March 15. doi: [10.1038/s41586-023-05795-1](https://doi.org/10.1038/s41586-023-05795-1)