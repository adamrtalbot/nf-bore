process BWA_INDEX {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'biocontainers/bwa:0.7.17--hed695b0_7' }"

    input:
    tuple val(meta), path(fasta), path(fai)

    output:
    tuple val(meta), path(fasta), path("*.{amb,ann,bwt,pac,sa}"), emit: index
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bwa index \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt
    touch ${fasta}.pac
    touch ${fasta}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}