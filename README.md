# CO-SPATIAL: A Nextflow pipeline to integrate spatial-ATAC-RNA seq data

Co-spatial is a Nextflow pipeline, inspired by [nf-core](https://nf-co.re) practices, to generate spatial
ATAC-seq and spatial RNA-seq-ready counts & peaks outputs for joint spatial ATAC-RNA sequencing workflows.

### Requirements
[Nextflow](https://www.nextflow.io/docs/stable/install.html)

[Docker](https://docs.docker.com/engine/install/)

### Running the pipeline
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

### Future updates
- Joint CUT&Tag-RNA sequencing
- Integrate downstream analysis for spatial multi-omics