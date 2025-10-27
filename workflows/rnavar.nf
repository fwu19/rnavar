/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// local
include { GTF2BED                   } from '../modules/local/gtf2bed'
include { FASTP                   } from '../modules/local2/fastp.nf'
include { WRITE_CSV                   } from '../modules/local2/write_csv.nf'
include { ADD_METADATA                   } from '../modules/local2/add_metadata.nf'

// nf-core
include { CAT_FASTQ                 } from '../modules/nf-core/cat/fastq'
include { FASTQC                    } from '../modules/nf-core/fastqc'
include { FASTQC as FASTQC_TRIM     } from '../modules/nf-core/fastqc'
include { GATK4_BASERECALIBRATOR    } from '../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_BEDTOINTERVALLIST   } from '../modules/nf-core/gatk4/bedtointervallist'
include { GATK4_INDEXFEATUREFILE    } from '../modules/nf-core/gatk4/indexfeaturefile'
include { GATK4_INTERVALLISTTOOLS   } from '../modules/nf-core/gatk4/intervallisttools'
include { MULTIQC                   } from '../modules/nf-core/multiqc'
include { SAMTOOLS_INDEX            } from '../modules/nf-core/samtools/index'
include { SEQ2HLA                   } from '../modules/nf-core/seq2hla/main'
include { TABIX_TABIX as TABIX      } from '../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIXGVCF  } from '../modules/nf-core/tabix/tabix'
include { UMITOOLS_EXTRACT          } from '../modules/nf-core/umitools/extract'

// local
include { RECALIBRATE               } from '../subworkflows/local/recalibrate'
include { SPLITNCIGAR               } from '../subworkflows/local/splitncigar'
include { PREPARE_ALIGNMENT         } from '../subworkflows/local/prepare_alignment'

// nf-core
include { BAM_MARKDUPLICATES_PICARD } from '../subworkflows/nf-core/bam_markduplicates_picard'
include { FASTQ_ALIGN_STAR          } from '../subworkflows/nf-core/fastq_align_star'

// local
include { checkSamplesAfterGrouping } from '../subworkflows/local/utils_nfcore_rnavar_pipeline'

/*
========================================================================================
    RUN MAIN WORKFLOW RNAVAR
========================================================================================
*/

workflow RNAVAR {
    take:
    input
    dbsnp
    dbsnp_tbi
    dict
    exon_bed
    fasta
    fasta_fai
    gtf
    known_sites
    known_sites_tbi
    star_index
    seq_center
    seq_platform
    aligner
    bam_csi_index
    extract_umi
    generate_gvcf
    skip_multiqc
    skip_baserecalibration
    skip_intervallisttools
    star_ignore_sjdbgtf
    tools

    main:

    // To gather all QC reports and versions for MultiQC
    reports = Channel.empty()
    versions = Channel.empty()

    // Parse the input data
    parsed_input = input
        .groupTuple()
        .map { samplesheet -> checkSamplesAfterGrouping(samplesheet) }
        .branch { meta, fastqs, bam, bai, cram, crai, vcf, tbi ->
            single: fastqs.size() == 1
            return [meta, fastqs.flatten()]
            multiple: fastqs.size() > 1
            return [meta, fastqs.flatten()]
            bam: bam
            return [meta, bam, bai]
            cram: cram
            return [meta, cram, crai]
            vcf: vcf
            return [meta, vcf, tbi]
        }

    // MODULE: Prepare the alignment files (index BAM/CRAM files that are missing an index)
    PREPARE_ALIGNMENT(
        parsed_input.cram,
        parsed_input.bam,
    )
    versions = versions.mix(PREPARE_ALIGNMENT.out.versions)

    // MODULE: Concatenate FastQ files from same sample if required
    CAT_FASTQ(parsed_input.multiple)

    def cat_fastq = CAT_FASTQ.out.reads.mix(parsed_input.single)

    versions = versions.mix(CAT_FASTQ.out.versions)

    // MODULE: Generate QC summary using FastQC
    FASTQC(cat_fastq)
    reports = reports.mix(FASTQC.out.zip.collect { _meta, logs -> logs })
    versions = versions.mix(FASTQC.out.versions)

    //MODULE: Trim adapter sequences
    // !! works for pair ends only!!
    FASTP(cat_fastq)
    def trim_fastq = FASTP.out.reads

    // MODULE: Generate QC summary using FastQC
    FASTQC_TRIM(trim_fastq)
    reports = reports.mix(FASTQC_TRIM.out.zip.collect { _meta, logs -> logs })
    versions = versions.mix(FASTQC_TRIM.out.versions)

    // MODULE: Extract UMIs from reads

    def umi_extracted_reads = Channel.empty()
    if (extract_umi) {
        UMITOOLS_EXTRACT(
            trim_fastq
        )
        versions = versions.mix(UMITOOLS_EXTRACT.out.versions)
        umi_extracted_reads = UMITOOLS_EXTRACT.out.reads
    }
    else {
        umi_extracted_reads = trim_fastq
    }

    // MODULE: Prepare the interval list from the GTF file using GATK4 BedToIntervalList

    GATK4_BEDTOINTERVALLIST(exon_bed, dict)
    def interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    versions = versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

    // MODULE: Scatter one interval-list into many interval-files using GATK4 IntervalListTools
    def interval_list_split = Channel.empty()
    if (!skip_intervallisttools) {
        GATK4_INTERVALLISTTOOLS(interval_list)
        interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map { _meta, bed -> [bed] }.collect()
        versions = versions.mix(GATK4_INTERVALLISTTOOLS.out.versions)
    }
    else {
        interval_list_split = interval_list.map { _meta, bed -> bed }
    }

    // MODULE: HLATyping with Seq2HLA
    if (tools.contains('seq2hla')) {
        SEQ2HLA(umi_extracted_reads)
        versions = versions.mix(SEQ2HLA.out.versions)
    }

    // SUBWORKFLOW: Perform read alignment using STAR aligner

    if (aligner == 'star') {
        FASTQ_ALIGN_STAR(
            umi_extracted_reads,
            star_index,
            gtf,
            star_ignore_sjdbgtf,
            seq_platform,
            seq_center,
            fasta,
            [[:], []],
        )
        //transcripts_fasta)

        def genome_bam = FASTQ_ALIGN_STAR.out.bam

        // Gather QC reports
        reports = reports.mix(FASTQ_ALIGN_STAR.out.log_out.collect { _meta, log_out -> log_out })
        reports = reports.mix(FASTQ_ALIGN_STAR.out.log_final.collect { it[1] }.ifEmpty([]))
        versions = versions.mix(FASTQ_ALIGN_STAR.out.versions)

        // SUBWORKFLOW: Mark duplicates with GATK4
        BAM_MARKDUPLICATES_PICARD(
            genome_bam,
            fasta,
            fasta_fai,
        )

        def markduplicate_indices = BAM_MARKDUPLICATES_PICARD.out.bai
            .mix(BAM_MARKDUPLICATES_PICARD.out.csi)
            .mix(BAM_MARKDUPLICATES_PICARD.out.crai)

        def genome_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam
            .join(markduplicate_indices, failOnDuplicate: true, failOnMismatch: true)
            .mix(PREPARE_ALIGNMENT.out.bam)

        //Gather QC reports
        reports = reports.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect { it[1] }.ifEmpty([]))
        reports = reports.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect { it[1] }.ifEmpty([]))
        reports = reports.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect { it[1] }.ifEmpty([]))
        reports = reports.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect { it[1] }.ifEmpty([]))
        versions = versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

        // SUBWORKFLOW: SplitNCigarReads from GATK4 over the intervals
        // Splits reads that contain Ns in their cigar string(e.g. spanning splicing events in RNAseq data).

        SPLITNCIGAR(
            genome_bam_bai,
            fasta,
            fasta_fai,
            dict,
            interval_list_split,
        )

        def splitncigar_bam_bai = SPLITNCIGAR.out.bam_bai
        versions = versions.mix(SPLITNCIGAR.out.versions)

        // MODULE: BaseRecalibrator from GATK4
        // Generates a recalibration table based on various co-variates
        def bam_variant_calling = Channel.empty()

        if (!skip_baserecalibration) {
            def interval_list_recalib = interval_list.map { _meta, bed -> [bed] }.flatten()
            def splitncigar_bam_bai_interval = splitncigar_bam_bai.combine(interval_list_recalib)

            GATK4_BASERECALIBRATOR(
                splitncigar_bam_bai_interval,
                fasta,
                fasta_fai,
                dict,
                known_sites,
                known_sites_tbi,
            )
            def bqsr_table = GATK4_BASERECALIBRATOR.out.table

            // Gather QC reports
            reports = reports.mix(bqsr_table.map { _meta, table -> table })
            versions = versions.mix(GATK4_BASERECALIBRATOR.out.versions)

            def bam_applybqsr = splitncigar_bam_bai.join(bqsr_table)

            def interval_list_applybqsr = interval_list.map { _meta, bed -> [bed] }.flatten()
            def applybqsr_bam_bai_interval = bam_applybqsr
                .combine(interval_list_applybqsr)
                .map { meta, bam, bai, table, interval -> [meta, bam, bai, table, interval] }

            // MODULE: ApplyBaseRecalibrator from GATK4
            // Recalibrates the base qualities of the input reads based on the recalibration table produced by the GATK BaseRecalibrator tool.
            RECALIBRATE(
                skip_multiqc,
                applybqsr_bam_bai_interval,
                dict.map { _meta, dict_ -> [dict_] },
                fasta_fai.map { _meta, fai -> fai },
                fasta.map { _meta, fasta_ -> [fasta_] },
            )

            bam_variant_calling = RECALIBRATE.out.bam

            // Gather QC reports
            reports = reports.mix(RECALIBRATE.out.qc.collect { it[1] }.ifEmpty([]))
            versions = versions.mix(RECALIBRATE.out.versions)
        }
        else {
            bam_variant_calling = splitncigar_bam_bai
        }


        // write out recalibrated.csv
        def my_dir = new File("${params.outdir}")
        def outdir = my_dir.absolutePath
        WRITE_CSV(
            bam_variant_calling
                .map { meta, bam, bai ->
                [ [ patient: meta.id ] + [ sample: meta.id ] + [ bam: "${outdir}/preprocessing/" + meta.id + "/" + bam.name ] + [ bai: "${outdir}/preprocessing/" + meta.id + "/" + bai.name ] ]
                }
                .collect(),
            "recalibrated.csv"
        )
        csv = WRITE_CSV.out.csv

        if (params.metadata){
            ADD_METADATA(
                csv,
                file(params.metadata, checkIfExists: true),
                'recalibrated.csv'
            )
        }

    }

    emit:
    reports  = reports // channel: qc reports for multiQC
    versions = versions // channel: [ path(versions.yml) ]
}
