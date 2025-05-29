process PLOT_VARIANT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-gviz=1.46.1 bioconda::bioconductor-biostrings=2.70.1 bioconda::bioconductor-genomicranges=1.54.1 bioconductor-variantannotation=1.48.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:b6524911af823c7c52518f6c886e5d1cf1f2cd6f-0' :
        'quay.io/biocontainers/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:b6524911af823c7c52518f6c886e5d1cf1f2cd6f-0' }"

    input:
    tuple val(meta), path(fasta), path(fai), path(bam), path(bam_index), path(gff), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plot_from = task.ext.plot_from ?: "27500"
    def plot_to = task.ext.plot_to ?: "28200"
    """
    #!/usr/bin/env Rscript --vanilla
    library(Gviz)
    library(GenomicRanges)
    library(Biostrings)
    library(VariantAnnotation)

    # Genome
    genome <- readDNAStringSet(file = "$fasta")
    gtrack <- GenomeAxisTrack(genome@ranges, name = "Genome")
    strack <- SequenceTrack(genome, name = "Sequence")

    # GFF
    gff <- read.table("$gff")
    gff <- gff[!gff\$V3 %in% c("gene",
                        "transcript",
                        "stop_codon",
                        "five_prime_UTR",
                        "three_prime_UTR"), ]
    gff <- gff[!gff\$V2 %in% c("ensembl_havana"), ]
    write.table(
        gff,
        file = "tmp.gff",
        quote = FALSE,
        col.names = FALSE,
        row.names = FALSE,
        sep = "\\t"
    )
    ftrack <-
        GeneRegionTrack(
            "tmp.gff",
            name = "Gene",
            showId = TRUE
        )

    # Alignment track
    atrack <-
        AlignmentsTrack(
            "$bam",
            isPaired = TRUE,
            name = "Alignments"
        )

    # Variants track
    vtrack <-
        AnnotationTrack(readVcf("$vcf")@rowRanges)

    # Plot
    pdf(file="${prefix}.pdf")
    plotTracks(
        list(strack, gtrack, ftrack, atrack, vtrack),
        from = $plot_from,
        to = $plot_to,
        showId = TRUE,
        add53 = TRUE,
        add35 = TRUE,
        labelPos = "below",
        stacking = "full"
    )
    dev.off()

    # Create versions file
    writeLines(c('"${task.process}":',
                 paste0('    R: "', paste0(R.Version()[c("major","minor")], collapse = "."), '"'),
                 paste0('    Gviz: "', packageVersion("Gviz"), '"'),
                 paste0('    GenomicRanges: "', packageVersion("GenomicRanges"), '"'),
                 paste0('    Biostrings: "', packageVersion("Biostrings"), '"'),
                 paste0('    VariantAnnotation: "', packageVersion("VariantAnnotation"), '"')),
               "versions.yml")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //; s/ (.*//')
    END_VERSIONS
    """
}