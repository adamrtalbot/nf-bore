#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-bore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Extract and analyze specific genomic regions from BAM files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/adamrtalbot/nf-bore
----------------------------------------------------------------------------------------
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NFBORE } from './workflows/nfbore'

//
// WORKFLOW: Run main nf-bore analysis pipeline
//
workflow NFCORE_NFBORE {
    NFBORE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    
    //
    // Input validation and channel creation
    //
    
    // Check mandatory parameters
    if (!params.fasta) {
        error "Please provide a reference FASTA file with --fasta"
    }
    if (!params.bams) {
        error "Please provide BAM files with --bams"
    }
    if (!params.gff) {
        error "Please provide a GFF annotation file with --gff"
    }
    if (!params.region) {
        error "Please provide genomic regions with --region"
    }

    //
    // Create input channels
    //
    
    // FASTA channel
    ch_fasta = Channel.fromFilePairs(params.fasta + "{,.fai}", checkIfExists: true, flat: true)
        .map { id, fasta, fai -> 
            def parsed_id = file(id).simpleName
            [ [fasta: parsed_id], fasta, fai ]
        }

    // BAM channel  
    ch_bam = Channel.fromFilePairs(params.bams + "{,.*ai}", flat: true, checkIfExists: true)
        .map { id, bam, bai -> 
            [ [id: id], bam, bai ]
        }

    // GFF channel
    ch_gff = Channel.fromPath(params.gff, checkIfExists: true)
        .map { gff ->
            [ [gff: gff.simpleName], gff ]
        }

    // Region channel
    ch_region = Channel.fromList(params.region.tokenize(','))

    //
    // WORKFLOW: Run pipeline
    //
    NFBORE (
        ch_fasta,
        ch_bam,
        ch_gff,
        ch_region
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/