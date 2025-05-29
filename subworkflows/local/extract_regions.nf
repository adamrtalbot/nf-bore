//
// Extract genomic regions from BAM, FASTA, and GFF files
//

include { EXTRACT_REGION_FROM_BAM   } from '../../modules/local/extract_region_from_bam'
include { EXTRACT_REGION_FROM_FASTA } from '../../modules/local/extract_region_from_fasta'
include { EXTRACT_REGION_FROM_GFF   } from '../../modules/local/extract_region_from_gff'
include { TRANSPOSE_GFF             } from '../../modules/local/transpose_gff'

workflow EXTRACT_REGIONS {

    take:
    ch_input // channel: [ val(meta), path(bam), path(bai), path(fasta), path(fai), path(gff) ]

    main:

    ch_versions = Channel.empty()

    //
    // Extract region from BAM
    //
    ch_bam_input = ch_input.map { meta, bam, bai, fasta, fai, gff ->
        [ meta, bam, bai, fasta, fai ]
    }
    EXTRACT_REGION_FROM_BAM ( ch_bam_input )
    ch_versions = ch_versions.mix(EXTRACT_REGION_FROM_BAM.out.versions.first())

    //
    // Extract region from FASTA
    //
    ch_fasta_input = ch_input.map { meta, bam, bai, fasta, fai, gff ->
        [ meta, fasta, fai ]
    }
    EXTRACT_REGION_FROM_FASTA ( ch_fasta_input )
    ch_versions = ch_versions.mix(EXTRACT_REGION_FROM_FASTA.out.versions.first())

    //
    // Extract region from GFF and transpose coordinates
    //
    ch_gff_input = ch_input.map { meta, bam, bai, fasta, fai, gff ->
        [ meta, gff ]
    }
    EXTRACT_REGION_FROM_GFF ( ch_gff_input )
    ch_versions = ch_versions.mix(EXTRACT_REGION_FROM_GFF.out.versions.first())

    TRANSPOSE_GFF ( EXTRACT_REGION_FROM_GFF.out.gff )
    ch_versions = ch_versions.mix(TRANSPOSE_GFF.out.versions.first())

    emit:
    bam     = EXTRACT_REGION_FROM_BAM.out.bam     // channel: [ val(meta), path(bam), path(bai) ]
    fasta   = EXTRACT_REGION_FROM_FASTA.out.fasta // channel: [ val(meta), path(fasta), path(fai) ]
    gff     = TRANSPOSE_GFF.out.gff               // channel: [ val(meta), path(gff) ]
    versions = ch_versions                        // channel: [ versions.yml ]
}