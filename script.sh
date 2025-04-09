rm output.out
rm output.err
module load nextflow-23.10.0
module load singularityce-4.0.3

bsub -n4 -R 'select[mem>4000] rusage[mem=4000]' -M4000 -G team353 -o output.out -e output.err \
    nextflow run main.nf \
    -profile singularity \
    --skip_fetch_data true \
    --input data/samplesheet.csv \
    --fastq_out data \
    --outdir result \
    -resume 2>&1 | tee headnode.log