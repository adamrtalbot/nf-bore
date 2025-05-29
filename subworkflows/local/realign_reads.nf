//
// Re-align reads to extracted genomic regions
//

include { EXTRACT_FASTQ_FROM_BAM } from '../../modules/local/extract_fastq_from_bam'
include { BWA_INDEX              } from '../../modules/local/bwa_index'
include { BWA_MEM                } from '../../modules/local/bwa_mem'

workflow REALIGN_READS {

    take:
    ch_region_bam   // channel: [ val(meta), path(bam), path(bai) ]
    ch_region_fasta // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    ch_versions = Channel.empty()

    //
    // Extract FASTQ from region BAM
    //
    EXTRACT_FASTQ_FROM_BAM ( ch_region_bam )
    ch_versions = ch_versions.mix(EXTRACT_FASTQ_FROM_BAM.out.versions.first())

    //
    // Create BWA index for region FASTA
    //
    BWA_INDEX ( ch_region_fasta )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    //
    // Combine FASTQ with FASTA and index for alignment
    //
    ch_align_input = EXTRACT_FASTQ_FROM_BAM.out.fastq
        .join(ch_region_fasta)
        .join(BWA_INDEX.out.index)
        .map { meta, fastq1, fastq2, fasta, fai, fasta_idx, index ->
            [ meta, [fastq1, fastq2], fasta_idx, index ]
        }

    //
    // Re-align reads to extracted region
    //
    BWA_MEM ( ch_align_input )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    emit:
    bam      = BWA_MEM.out.bam      // channel: [ val(meta), path(bam), path(bai) ]
    versions = ch_versions          // channel: [ versions.yml ]
}