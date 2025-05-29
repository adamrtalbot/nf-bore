process TRANSPOSE_GFF {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bioconductor-gviz=1.46.1 bioconda::bioconductor-biostrings=2.70.1 bioconda::bioconductor-genomicranges=1.54.1 bioconductor-variantannotation=1.48.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:b6524911af823c7c52518f6c886e5d1cf1f2cd6f-0' :
        'quay.io/biocontainers/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:b6524911af823c7c52518f6c886e5d1cf1f2cd6f-0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.transposed.gff3"), emit: gff
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr = meta.region.split(":")[0]
    def start = meta.region.split(":")[1].split("-")[0]
    def stop = meta.region.split(":")[1].split("-")[1]
    """
    #!/usr/bin/env Rscript
    gff <- read.table("$gff")
    gff\$V1 <- "${meta.region}"
    gff\$V4 <- gff\$V4 - $start
    gff\$V5 <- gff\$V5 - $start
    gff <- gff[!gff\$V3 %in% c("gene",
                        "transcript",
                        "stop_codon",
                        "five_prime_UTR",
                        "three_prime_UTR"), ]
    gff <- gff[!gff\$V2 %in% c("ensembl_havana"), ]
    write.table(gff, 
                file = "${prefix}.transposed.gff3", 
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\\t")

    # Create versions file
    writeLines(c('"${task.process}":',
                 paste0('    R: "', paste0(R.Version()[c("major","minor")], collapse = "."), '"')),
               "versions.yml")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transposed.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //; s/ (.*//')
    END_VERSIONS
    """
}