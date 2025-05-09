nextflow run main.nf \
    -profile docker \
    --skip_fetch_data true \
    --skip_star_genome true \
    --input data/samplesheet.csv \
    --fastq_out data \
    --spatial_barcodes /mnt/raid5/kbui/misc/co-spatial/assets/spatial_barcodes.txt \
    --ref_map data/references/Mus_musculus/Ensembl/GRCm38/Sequence/star \
    --ref_annotation data/references/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf \
    --ref_atac_genome /mnt/raid5/kbui/misc/co-spatial/data/references/Mus_musculus/mm10-2020-A/mm10 \
    --outdir result \
    -resume 2>&1 | tee headnode.log