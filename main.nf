process EXTRACT_REGION_FROM_BAM {

    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    conda "bioconda::samtools=1.19.2"

    input:
        tuple val(meta), path(bam), path(index), path(fasta), path(fai) 

    output:
        tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam

    script:
    """
    samtools view \\
        -O BAM \\
        -T $fasta \\
        $bam \\
        ${meta.region} \\
    | samtools sort > ${meta.prefix}.bam \\
    && samtools index ${meta.prefix}.bam
    """
}

process EXTRACT_REGION_FROM_FASTA {

    container "quay.io/biocontainers/mulled-v2-f01e242bdea19948f0576fdca94777242fe4c2cb:4238fb992d2a93e648108c86e3a9f51348e834a9-0" // Mulled with samtools and bedtools
    conda "bioconda::samtools=1.19.2 bioconda::bedtools=2.31.1"

    input:
        tuple val(meta), path(fasta), path(fai)

    output:
        tuple val(meta), path("*.region.fasta"), path("*.region.fasta.fai"), emit: fasta

    script:
    // TODO: Sort this out
    def region_tab = meta.region.replaceAll(":", '\t').replaceAll("-", '\t') + "\t" + meta.region.replaceAll(":", '_').replaceAll("-", '_')
    """
    bedtools getfasta -nameOnly -fi $fasta -bed <(echo "$region_tab") > ${meta.prefix}.region.fasta && samtools faidx ${meta.prefix}.region.fasta
    """
}

process EXTRACT_REGION_FROM_GFF {

    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"
    conda "bioconda::bedtools=2.31.1"

    input:
        tuple val(meta), path(gff)

    output:
        tuple val(meta), path("*.gff3"), emit: gff

    script:
    def region_tab = meta.region.replaceAll(":", '\t').replaceAll("-", '\t')
    """
    bedtools intersect -a $gff -b <(echo "$region_tab") > ${meta.prefix}.gff3
    """
}

process TRANSPOSE_GFF {

    // tmp Wave container. Will expire.
    container "195996028523.dkr.ecr.eu-west-1.amazonaws.com/wave/build:bioconductor-gviz-1.46.1_bioconductor-biostrings-2.70.1_bioconductor-genomicrang--c0e7a90fe5f0eba3"
    // Lets put this in an environment.yml...
    conda "bioconda::bioconductor-gviz=1.46.1 bioconda::bioconductor-biostrings=2.70.1 bioconda::bioconductor-genomicranges=1.54.1 bioconductor-variantannotation=1.48.1"

    input:
        tuple val(meta), path(gff)

    output:
        tuple val(meta), path("*.gff3"), emit: gff

    script:
    def chr   = meta.region.split(":")[0] // Matches output of bedtools intersect
    def start = meta.region.split(":")[1].split("-")[0]
    def stop  = meta.region.split(":")[1].split("-")[1]
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
    write.table(gff, file = "${meta.prefix}.transposed.gff3", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\\t")
    """
}

process EXTRACT_FASTQ_FROM_BAM {

    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    conda "bioconda::samtools=1.19.2"

    input:
        tuple val(meta), path(bam), path(index)
    
    output:
        tuple val(meta), path("*.1.fastq.gz"), path("*.2.fastq.gz"), emit: fastq

    script:
    """
    samtools collate -u -O $bam | samtools fastq -1 ${meta.prefix}.1.fastq.gz -2 ${meta.prefix}.2.fastq.gz -0 /dev/null -s /dev/null -n
    """
}

process CREATE_BWA_INDEX {

    container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    conda "bioconda::bwa=0.7.17"

    input:
        tuple val(meta), path(fasta), path(fai)

    output:
        tuple val(meta), path("${fasta_name}.bwt"), path("${fasta_name}.pac"), path("${fasta_name}.amb"), path("${fasta_name}.sa"), path("${fasta_name}.ann"), emit: bwa_index

    script:
    fasta_name = fasta.name
    """
    bwa index -p $fasta_name $fasta
    """
}

process ALIGN_WITH_BWA {

    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:a34558545ae1413d94bde4578787ebef08027945-0"
    conda "bioconda::bwa=0.7.17 bioconda::samtools=1.19.2"

    input:
        tuple val(meta), path(fastq1), path(fastq2), path(fasta), path(fai), path(bwt), path(pac), path(amb), path(sa), path(ann)

    output:
        tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam

    script:
    """
    bwa mem \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}" \\
        $fasta \\
        $fastq1 \\
        $fastq2 \\
    | samtools sort \\
        -O BAM \\
        > ${meta.prefix}.bam \\
    && samtools index ${meta.prefix}.bam
    """
}

process CREATE_SEQUENCE_DICTIONARY {
    
    container "quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0"
    conda "bioconda::gatk4=4.5.0.0"

    input:
        tuple val(meta), path(fasta), path(fai)

    output:
        tuple val(meta), path("*.dict"), emit: fasta_dict

    script:
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[GATK CreateSequenceDictionary] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CreateSequenceDictionary \\
        --REFERENCE $fasta \\
        --URI $fasta \\
        --TMP_DIR .
    """
}

process CALL_HAPLOTYPES {

    container "quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0"
    conda "bioconda::gatk4=4.5.0.0"

    input:
        tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(dict)

    output:
        tuple val(meta), path("*.vcf.gz"), emit: vcf

    script:
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        HaplotypeCaller \\
        --input $bam \\
        --output ${meta.prefix}.vcf.gz \\
        --reference $fasta \\
        --tmp-dir .
    """
}

process INDEX_VCF {
    
        container "quay.io/biocontainers/htslib:1.19.1--h81da01d_1"
        conda "bioconda::htslib=1.19.1"
    
        input:
            tuple val(meta), path(vcf)
    
        output:
            tuple val(meta), path("*.tbi"), emit: vcf_index
    
        script:
        """
        tabix $vcf
        """
}

process PLOT_VARIANT {

    container ""
    conda "bioconda::bioconductor-gviz=1.46.1 bioconda::bioconductor-biostrings=2.70.1 bioconda::bioconductor-genomicranges=1.54.1 bioconductor-variantannotation=1.48.1"

    cpus 2
    memory '4 GB'

    input:
        tuple val(meta), path(fasta), path(fai), path(bam), path(bam_index), path(gff), path(vcf), path(vcf_index)

    output:
        path "*.pdf"

    script:
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
    pdf(file="${meta.prefix}.pdf")
    plotTracks(
        list(strack, gtrack, ftrack, atrack, vtrack),
        from = 27500,
        to = 28200,
        showId = TRUE,
        add53 = TRUE,
        add35 = TRUE,
        labelPos = "below",
        stacking = "full"
    )
    dev.off()
    """
}

workflow {
    ch_fasta = Channel.fromFilePairs(params.fasta + "{,.fai}", checkIfExists: true, flat: true)
                    .map { id, fasta, fai -> 
                        def parsed_id = file(id).simpleName
                        tuple([fasta: parsed_id], fasta, fai)
                    }
    ch_fasta.view()
    ch_bam   = Channel.fromFilePairs(params.bams + "{,.*ai}", flat: true  , checkIfExists: true)
                    .map { id, bam, bai -> 
                        tuple([id: id], bam, bai)
                    }
    ch_bam.view()
    ch_gff   = Channel.fromPath(params.gff, checkIfExists: true)
                    .map { gff ->
                        tuple([gff: gff.simpleName], gff)
                    }
    ch_gff.view()
    ch_region = Channel.fromList(params.region.tokenize(','))
    ch_region.view()

    // First, combine all inputs into a single big input channel containing all combinations
    ch_input = ch_bam
        .combine(ch_fasta)
        .map { bam_meta, bam, bai, fasta_meta, fasta, fai -> 
             tuple(bam_meta + fasta_meta, bam, bai, fasta, fai)
        }
        .combine(ch_gff)
        .map { meta, bam, bai, fasta, fai, gff_meta, gff ->
             tuple(meta + gff_meta, bam, bai, fasta, fai, gff)
        }
        .combine ( ch_region )
        .map { meta, bam, bai, fasta, fai, gff, region ->
             tuple(meta + [region: region], bam, bai, fasta, fai, gff)
        }
        // Add a prefix to use for file names
        .map { meta, bam, bai, fasta, fai, gff->
            tuple(
                meta + [ prefix: [
                    meta.id, 
                    meta.region.replaceAll(":", '_').replaceAll("-", '_')
                ].join('_') ],
                bam, 
                bai, 
                fasta, 
                fai, 
                gff
            )
        }


    region_bam = ch_input
        .map { meta, bam, bai, fasta, fai, gff ->
            tuple(meta, bam, bai, fasta, fai)
        }
        | EXTRACT_REGION_FROM_BAM

    region_fastq = region_bam
        | EXTRACT_FASTQ_FROM_BAM

    region_fasta = ch_input
        .map { meta, bam, bai, fasta, fai, gff ->
            tuple(meta, fasta, fai)
        }
        | EXTRACT_REGION_FROM_FASTA

    region_index = region_fasta
        | CREATE_BWA_INDEX

    realigned_bam = region_fastq
        .join(region_fasta)
        .join(region_index)
        | ALIGN_WITH_BWA

    region_gff = ch_input
        .map { meta, bam, bai, fasta, fai, gff ->
            tuple(meta, gff)
        }
        | EXTRACT_REGION_FROM_GFF
        | TRANSPOSE_GFF

    region_dict = region_fasta
        | CREATE_SEQUENCE_DICTIONARY

    region_variants = realigned_bam
        .join(region_fasta)
        .join(region_dict)
        | CALL_HAPLOTYPES

    region_variants_index = INDEX_VCF(region_variants)

    ch_input
        .map { it[0] }
        .join(region_fasta)
        .join(realigned_bam)
        .join(region_gff)
        .join(region_variants)
        .join(region_variants_index)
        | PLOT_VARIANT

}