nextflow.enable.dsl2 = 2

process calculate_scrublet_score {
    conda ''

    input:

    output:

    script
    """
    """
}

workflow {
    Channel.fromPath(params.sample_file) \
    | splitCsv(header : true) \
    | map { row -> tuple(row.sample, row.directory)} \
    | calculate_scrublet_score 
}