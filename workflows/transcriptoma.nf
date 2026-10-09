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

    SAMTOOLS_SORT(HISAT2_ALIGN.out.bam)


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
