process EXTRACT_FASTQ_FROM_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.1.fastq.gz"), path("*.2.fastq.gz"), emit: fastq
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools collate \\
        -u \\
        -O \\
        $args \\
        $bam \\
    | samtools fastq \\
        -1 ${prefix}.1.fastq.gz \\
        -2 ${prefix}.2.fastq.gz \\
        -0 /dev/null \\
        -s /dev/null \\
        -n

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.fastq.gz
    touch ${prefix}.2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}