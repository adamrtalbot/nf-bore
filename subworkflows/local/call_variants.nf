//
// Call variants using GATK HaplotypeCaller
//

include { GATK_CREATESEQUENCEDICTIONARY } from '../../modules/local/gatk_createsequencedictionary'
include { GATK_HAPLOTYPECALLER          } from '../../modules/local/gatk_haplotypecaller'
include { TABIX_TABIX                   } from '../../modules/local/tabix_tabix'

workflow CALL_VARIANTS {

    take:
    ch_bam   // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    ch_versions = Channel.empty()

    //
    // Create sequence dictionary
    //
    GATK_CREATESEQUENCEDICTIONARY ( ch_fasta )
    ch_versions = ch_versions.mix(GATK_CREATESEQUENCEDICTIONARY.out.versions.first())

    //
    // Combine BAM with FASTA and dictionary for variant calling
    //
    ch_haplotypecaller_input = ch_bam
        .join(ch_fasta)
        .join(GATK_CREATESEQUENCEDICTIONARY.out.dict)
        .map { meta, bam, bai, fasta, fai, dict ->
            [ meta, bam, bai, fasta, fai, dict ]
        }

    //
    // Call variants
    //
    GATK_HAPLOTYPECALLER ( ch_haplotypecaller_input )
    ch_versions = ch_versions.mix(GATK_HAPLOTYPECALLER.out.versions.first())

    //
    // Index VCF
    //
    TABIX_TABIX ( GATK_HAPLOTYPECALLER.out.vcf )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    emit:
    vcf      = GATK_HAPLOTYPECALLER.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi      = TABIX_TABIX.out.tbi          // channel: [ val(meta), path(tbi) ]
    versions = ch_versions                  // channel: [ versions.yml ]
}