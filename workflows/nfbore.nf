//
// nf-bore: Extract and analyze specific genomic regions
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { EXTRACT_REGIONS } from '../subworkflows/local/extract_regions'
include { REALIGN_READS   } from '../subworkflows/local/realign_reads'
include { CALL_VARIANTS   } from '../subworkflows/local/call_variants'
include { PLOT_VARIANT    } from '../modules/local/plot_variant'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFBORE {

    take:
    ch_fasta  // channel: [ val(meta), path(fasta), path(fai) ]
    ch_bam    // channel: [ val(meta), path(bam), path(bai) ]
    ch_gff    // channel: [ val(meta), path(gff) ]
    ch_region // channel: [ region ]

    main:

    ch_versions = Channel.empty()

    //
    // Combine all inputs into a single channel with all combinations
    //
    ch_input = ch_bam
        .combine(ch_fasta)
        .map { bam_meta, bam, bai, fasta_meta, fasta, fai -> 
             [ bam_meta + fasta_meta, bam, bai, fasta, fai ]
        }
        .combine(ch_gff)
        .map { meta, bam, bai, fasta, fai, gff_meta, gff ->
             [ meta + gff_meta, bam, bai, fasta, fai, gff ]
        }
        .combine(ch_region)
        .map { meta, bam, bai, fasta, fai, gff, region ->
             [ meta + [region: region], bam, bai, fasta, fai, gff ]
        }
        // Add a prefix to use for file names
        .map { meta, bam, bai, fasta, fai, gff ->
            [ meta + [ prefix: [
                meta.id, 
                meta.region.replaceAll(":", '_').replaceAll("-", '_')
            ].join('_') ], bam, bai, fasta, fai, gff ]
        }

    //
    // SUBWORKFLOW: Extract regions from BAM, FASTA, and GFF
    //
    EXTRACT_REGIONS ( ch_input )
    ch_versions = ch_versions.mix(EXTRACT_REGIONS.out.versions)

    //
    // SUBWORKFLOW: Re-align reads to extracted regions
    //
    REALIGN_READS ( 
        EXTRACT_REGIONS.out.bam,
        EXTRACT_REGIONS.out.fasta
    )
    ch_versions = ch_versions.mix(REALIGN_READS.out.versions)

    //
    // SUBWORKFLOW: Call variants on re-aligned reads
    //
    CALL_VARIANTS (
        REALIGN_READS.out.bam,
        EXTRACT_REGIONS.out.fasta
    )
    ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)

    //
    // MODULE: Create visualization plots
    //
    ch_plot_input = ch_input
        .map { meta, bam, bai, fasta, fai, gff -> meta }
        .join(EXTRACT_REGIONS.out.fasta)
        .join(REALIGN_READS.out.bam)
        .join(EXTRACT_REGIONS.out.gff)
        .join(CALL_VARIANTS.out.vcf)
        .join(CALL_VARIANTS.out.tbi)
        .map { meta, fasta, fai, bam, bai, gff, vcf, tbi ->
            [ meta, fasta, fai, bam, bai, gff, vcf, tbi ]
        }

    PLOT_VARIANT ( ch_plot_input )
    ch_versions = ch_versions.mix(PLOT_VARIANT.out.versions)

    emit:
    region_bam    = EXTRACT_REGIONS.out.bam    // channel: [ val(meta), path(bam), path(bai) ]
    region_fasta  = EXTRACT_REGIONS.out.fasta  // channel: [ val(meta), path(fasta), path(fai) ]
    region_gff    = EXTRACT_REGIONS.out.gff    // channel: [ val(meta), path(gff) ]
    realigned_bam = REALIGN_READS.out.bam      // channel: [ val(meta), path(bam), path(bai) ]
    vcf           = CALL_VARIANTS.out.vcf      // channel: [ val(meta), path(vcf) ]
    plots         = PLOT_VARIANT.out.pdf       // channel: [ val(meta), path(pdf) ]
    versions      = ch_versions                // channel: [ versions.yml ]
}