params.input_bcf   = '/home/zack/experiment_6_3_1/bcf_filter/wgs_cohort.bcf'
params.chunks_file = 'chunks_list.txt'
params.min_qual    = 30
params.min_depth   = 10
params.chunk_size  = 1000

process filterVariants {
    publishDir '/home/zack/experiment_6_3_1/bcf_filter/', mode: 'copy', overwrite: true

    input:
    tuple path(bcf_file), val(chunk)

    output:
    path "filtered_${chunk.replace(':','_').replace('-','_')}.bcf"

    script:
    """
    /home/zack/experiment_6_3_1/bcf_filter/target/release/bcf_filter_tool \\
        --input ${bcf_file} \\
        --output filtered_${chunk.replace(':','_').replace('-','_')}.bcf \\
        --min-qual ${params.min_qual} \\
        --min-depth ${params.min_depth} \\
        --chunk-size ${params.chunk_size} \\
        || { echo "Error processing chunk ${chunk}"; exit 1; }
    """
}

workflow {
    bcf_file_channel = Channel.fromPath(params.input_bcf)
    chunk_channel = Channel.from(['chr1:1-2000', 'chr2:1-3000'])
    filter_inputs = bcf_file_channel.combine(chunk_channel)
    filter_inputs.map { bcf_file, chunk -> tuple(bcf_file, chunk) }.set { filter_inputs_unpacked }
    filterVariants(filter_inputs_unpacked)
    filterVariants.out.view()
}