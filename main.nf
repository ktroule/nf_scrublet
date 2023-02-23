#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process calculate_scrublet_score {
    conda '/nfs/team292/kt22/pediatric_gonads/env/p-PG-scanpy1.9.1_20220517.txt'

    input:
    tuple val(sample), val(directory)

    script
    """
    00_calculate_scrublet_score.py \
        --sample $sample \
        --directory $directory \
        --min_genes 100 \
        --min_cells 100 \
        --rnd_seed 42 \
        --out_dir /nfs/team292/kt22/nf_tests
    """
}

workflow {
    Channel.fromPath(params.sample_file) \
        | splitCsv(header : true) \
        | map { row -> tuple(row.sample, row.directory) } \
        | calculate_scrublet_score 
}