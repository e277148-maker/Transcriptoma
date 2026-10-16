/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_transcriptoma_pipeline'

include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD } from '../modules/nf-core/hisat2/build/main'
include { HISAT2_ALIGN } from '../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { STRINGTIE_STRINGTIE } from '../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE } from '../modules/nf-core/stringtie/merge/main'
include { GFFREAD } from '../modules/nf-core/gffread/main'
include { KALLISTO_INDEX } from '../modules/nf-core/kallisto/index/main'
include { KALLISTO_QUANT } from '../modules/nf-core/kallisto/quant/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TRANSCRIPTOMA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    outdir

    main:

    def ch_versions = channel.empty()

    ch_reference_fasta = channel.fromPath(params.reference_fasta)

    ch_reference_gtf   = channel.fromPath(params.reference_gtf)

// Index the reference genome FASTA for SAMtools.
    ch_reference_fasta_with_meta = ch_reference_fasta
        .map { fasta -> tuple([id: 'reference'], fasta, []) }

    SAMTOOLS_FAIDX(ch_reference_fasta_with_meta, false)

    ch_reference_fai = SAMTOOLS_FAIDX.out.fai

    ch_reference_gtf_with_meta = ch_reference_gtf.map { gtf ->tuple([id: 'reference'], gtf)}

    HISAT2_EXTRACTSPLICESITES(ch_reference_gtf_with_meta)

    ch_reference = ch_reference_fasta
    .combine(ch_reference_gtf)
    .map { fasta, gtf ->
        tuple([id: 'reference'], fasta, gtf)
    }

    ch_reference_for_build = ch_reference
    .combine(HISAT2_EXTRACTSPLICESITES.out.txt, by: 0)
    .map { meta, fasta, gtf, splicesites ->
        tuple(meta, fasta, gtf, splicesites)
    }

    HISAT2_BUILD(ch_reference_for_build, '1 GB')

    ch_hisat2_index = HISAT2_BUILD.out.index.first()
    ch_splicesites  = HISAT2_EXTRACTSPLICESITES.out.txt.first()

    HISAT2_ALIGN(ch_samplesheet, ch_hisat2_index, ch_splicesites, false)

// Prepare the reference FASTA and its index for SAMtools.
    ch_reference_for_sort = ch_reference_fasta
        .combine(ch_reference_fai)
        .map { fasta, fai ->
            tuple([id: 'reference'], fasta, fai)
        }

// Sort the alignment files.
    SAMTOOLS_SORT(
        HISAT2_ALIGN.out.bam,
        ch_reference_for_sort,
        ''
    )

// Prepare sorted BAMs for StringTie.
// These are short-read data, so no long-read BAM is provided.
    ch_stringtie_input = SAMTOOLS_SORT.out.bam.map { meta, bam -> tuple(meta, bam, [])}

// Assemble transcripts using the reference annotation.
    STRINGTIE_STRINGTIE(ch_stringtie_input, '',ch_reference_gtf.first())

// Merge transcript assemblies from all samples.
    ch_stringtie_gtfs = STRINGTIE_STRINGTIE.out.transcript_gtf
        .map { meta, gtf -> tuple(meta, gtf) }
        .collect()
        .map { gtfs -> tuple([id: 'merged'], gtfs) }

    ch_reference_gtf_for_merge = ch_reference_gtf
        .map { gtf -> tuple([id: 'reference'], gtf) }

    STRINGTIE_MERGE(ch_stringtie_gtfs, ch_reference_gtf_for_merge)

// Extract transcript sequences from the merged annotation.
    ch_transcripts_fasta_input = STRINGTIE_MERGE.out.merged_gtf
        .map { meta, gtf -> tuple(meta, gtf) }

    GFFREAD(ch_transcripts_fasta_input, ch_reference_fasta.first())

// Build the Kallisto transcriptome index.
    ch_transcripts_fasta = GFFREAD.out.gffread_fasta

    KALLISTO_INDEX(ch_transcripts_fasta)

// Quantify transcript abundance for each sample.
    ch_kallisto_index_for_quant = KALLISTO_INDEX.out.index
        .combine( channel.value([]) )
        .map { meta, index, unused -> tuple(meta, index, [], []) }

    KALLISTO_QUANT(
        ch_samplesheet,
        ch_kallisto_index_for_quant,
        250,
        100
    )


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'transcriptoma_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
