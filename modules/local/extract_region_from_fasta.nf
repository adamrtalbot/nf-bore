process EXTRACT_REGION_FROM_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.19.2 bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f01e242bdea19948f0576fdca94777242fe4c2cb:4238fb992d2a93e648108c86e3a9f51348e834a9-0' :
        'quay.io/biocontainers/mulled-v2-f01e242bdea19948f0576fdca94777242fe4c2cb:4238fb992d2a93e648108c86e3a9f51348e834a9-0' }"

    input:
    tuple val(meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.region.fasta"), path("*.region.fasta.fai"), emit: fasta
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def region_tab = meta.region.replaceAll(":", '\t').replaceAll("-", '\t') + "\t" + meta.region.replaceAll(":", '_').replaceAll("-", '_')
    """
    bedtools getfasta \\
        -nameOnly \\
        -fi $fasta \\
        -bed <(echo "$region_tab") \\
        $args \\
        > ${prefix}.region.fasta

    samtools faidx ${prefix}.region.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.region.fasta
    touch ${prefix}.region.fasta.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}