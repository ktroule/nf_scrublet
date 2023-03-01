#!/usr/bin/env nextflow2
nextflow.enable.dsl=2

process calculate_scrublet_score {
    conda '/nfs/team292/kt22/pediatric_gonads/env/p-PG-scanpy1.9.1_20220517.txt'
    input:
        tuple val(sample), val(directory)
    output:
        file "${sample}.h5ad"
    script:
    """
    00_calculate_scrublet_score.py \
        --sample $sample \
        --directory $directory \
        --min_genes $params.min_genes \
        --min_cells $params.min_cells \
        --rnd_seed $params.rnd_seed \
        --out_dir $params.output_dir
    """
}

workflow {
    out_test = Channel.fromPath(params.sample_file) \
        | splitCsv(header : true) \
        | map { row -> tuple(row.sample, row.directory) } \
        | calculate_scrublet_score
}