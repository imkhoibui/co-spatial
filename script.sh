nextflow run main.nf \
    -profile conda \
    --input data/samplesheet.csv \
    --fastq_out data \
    --outdir result \
    -resume 2>&1 | tee headnode.log