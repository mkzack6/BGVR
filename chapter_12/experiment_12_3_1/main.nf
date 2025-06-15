nextflow.enable.dsl=2

params.sample_list = "samples.txt"
params.output_dir = "results"
params.region = "chr1:1-32"
params.mock = true

workflow {
    // Define the 'samples_ch' channel inside the workflow block
    samples_ch = Channel
        .fromPath(params.sample_list)
        .splitText()
        .map { it.trim() }

    // Pass the samples to the alignmentOrFetch process
    bam_ch = alignmentOrFetch(samples_ch)

    // Variant calling process
    vcf_ch = variantCalling(bam_ch)

    // Merge VCF files
    mergeVcfs(vcf_ch.collect())
}

process alignmentOrFetch {
    input:
    val sample_id

    output:
    path "${sample_id}.bam"

    script:
    """
    # Create a mock BAM file for testing
    dd if=/dev/urandom of=${sample_id}.bam bs=1M count=10
    """
}

process variantCalling {
    input:
    path bam

    output:
    path "variants_${bam.simpleName}.vcf"

    script:
    if (params.mock) {
        """
        # Mock variant calling - create a dummy VCF file for testing
        echo '##fileformat=VCFv4.2' > variants_${bam.simpleName}.vcf
        echo '##source=MockVariantCaller' >> variants_${bam.simpleName}.vcf
        echo '##contig=<ID=chr1,length=248956422>' >> variants_${bam.simpleName}.vcf
        echo '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${bam.simpleName}' >> variants_${bam.simpleName}.vcf
        echo 'chr1\t14653\t.\tA\tG\t100\tPASS\t.\tGT\t0/1' >> variants_${bam.simpleName}.vcf
        echo 'chr1\t14907\t.\tA\tG\t100\tPASS\t.\tGT\t0/1' >> variants_${bam.simpleName}.vcf
        echo 'chr1\t15211\t.\tG\tA\t100\tPASS\t.\tGT\t1/1' >> variants_${bam.simpleName}.vcf
        """
    } else {
        """
        # Run our Rust BAM counter tool
        ./bam_counter --bam ${bam} --region chr1:1-1000000
        """
    }
}

process mergeVcfs {
    publishDir params.output_dir, mode: 'copy'

    input:
    path vcfs

    output:
    path "merged.vcf"

    script:
    """
    # Mock merge - concatenate the VCF files
    cat ${vcfs[0]} > header.tmp
    grep "^#" header.tmp > merged.vcf
    for vcf in ${vcfs}; do
        grep -v "^#" \$vcf >> merged.vcf || true
    done
    """
} 