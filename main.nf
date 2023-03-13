#!/usr/bin/env nextflow2
nextflow.enable.dsl=2

process calculate_scrublet_score {
    publishDir "${params.output_dir}", mode: 'copy'
    conda '/nfs/team292/kt22/nf_scrublet/env/nf_scrublet.yml'
    input:
        tuple val(sample), val(directory)
    output:
        path "${sample}.h5ad"
    script:
    """
    00_calculate_scrublet_score.py \
        --sample "${sample}" \
        --directory "${directory}" \
        --min_genes ${params.min_genes} \
        --min_cells ${params.min_cells} \
        --rnd_seed ${params.rnd_seed} \
        --out_dir "${params.output_dir}"
    """
}

workflow {
    Channel.fromPath(params.sample_file) \
        | splitCsv(header : true) \
        | map { row -> tuple(row.sample, row.directory) } \
        | calculate_scrublet_score
}
// /software/team292/kt22/nextflow/nextflow run main.nf -with-conda true -params-file params.yaml -c nextflow.config