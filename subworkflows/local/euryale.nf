//
// LOCAL MODULES
//
include { FASTX_COLLAPSER } from '../../modules/local/fastx_toolkit/collapser'

//
// SUBWORKFLOWS
//

include { PREPROCESS } from './preprocess'
include { TAXONOMY } from './taxonomy'
include { ASSEMBLY } from './assembly'
include { ALIGNMENT } from './alignment'
include { FUNCTIONAL } from './functional'

//
// MODULES: Installed directly from nf-core/modules
//

include { GUNZIP                      } from '../../modules/nf-core/gunzip/main'
include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/nf-core/custom/dumpsoftwareversions/main'

workflow EURYALE {

    take:
    unmapped_reads

    main:
    if (params.reference_fasta == null && params.diamond_db == null) { exit 1, 'A reference fasta (--reference_fasta) or a DIAMOND db (--diamond_db) must be specified' }

    ch_versions = Channel.empty()
    ch_kaiju_db = Channel.value([ [id: "kaiju_db"], file(params.kaiju_db)])
    ch_reference_fasta = params.reference_fasta ? file(params.reference_fasta) : []
    ch_diamond_db = params.diamond_db ? file(params.diamond_db) : []
    ch_id_mapping = params.id_mapping ? file(params.id_mapping) : []

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    PREPROCESS (
        unmapped_reads
    )
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)

    PREPROCESS.out.merged_reads
        .set{ clean_reads }

    if (params.assembly_based) {
        ASSEMBLY (
            unmapped_reads
        )
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
    }

    GUNZIP (
        clean_reads
    )

    GUNZIP.out.gunzip
        .set { decompressed_reads }

    FASTX_COLLAPSER (
        decompressed_reads
    )

    ALIGNMENT (
        FASTX_COLLAPSER.out.collapsed,
        ch_reference_fasta,
        ch_diamond_db
    )
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    TAXONOMY (
        clean_reads,
        ch_kaiju_db
    )
    ch_versions = ch_versions.mix(TAXONOMY.out.versions)

    FUNCTIONAL (
        ALIGNMENT.out.alignments,
        ch_id_mapping
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //

    ch_multiqc_files  = Channel.empty()

    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_multiqc_files  = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files  = ch_multiqc_files.mix(PREPROCESS.out.multiqc_files.collect())
    ch_multiqc_files  = ch_multiqc_files.mix(TAXONOMY.out.kaiju_report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files  = ch_multiqc_files.mix(ALIGNMENT.out.multiqc_files.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}
