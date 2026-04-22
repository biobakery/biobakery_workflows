#!/usr/bin/env nextflow

/*
 * bioBakery Visualization Workflow — Nextflow DSL2
 *
 * Equivalent to the AnADAMA2 vis.py workflow.
 * Supports:
 *   - Local execution (with or without Docker / Singularity / Conda)
 *   - HPC cluster (SLURM, PBS/Torque, LSF, SGE)
 *   - Cloud (AWS Batch, Google Cloud Batch)
 *
 * Usage:
 *   nextflow run main.nf --input /path/to/workflow/output [options]
 *
 * For multiple projects in parallel (glob pattern):
 *   nextflow run main.nf --input "/data/*/wmgx_output" [options]
 *
 * Copyright (c) 2021 Harvard School of Public Health — MIT License
 */

nextflow.enable.dsl = 2

// ============================================================================
// Parameter defaults — mirrors vis.py defaults exactly
// ============================================================================
params.input                  = null       // required
params.output                 = "vis_output"
params.format                 = "pdf"      // "pdf" | "html"
params.project_name           = ""
params.author_name            = ""
params.header_image           = ""
params.introduction_text      = ""
params.use_template           = ""

params.input_metadata         = ""
params.input_file_type        = []
params.input_picard           = ""
params.input_picard_extension = "quality_by_cycle_metrics"

params.metadata_categorical   = []
params.metadata_continuous    = []
params.metadata_exclude       = []
params.max_missing            = "20.0"

params.min_abundance          = 0.01
params.min_samples            = 3

params.max_sets_heatmap       = 25
params.max_sets_barplot       = 15
params.max_groups_barplot     = 5
params.correlation_threshold  = 0.7

// ============================================================================
// Validate required parameters
// ============================================================================
if (!params.input) {
    log.error """
    ═══════════════════════════════════════════════════════════════
    ERROR: --input is required.

    Provide the path to a biobakery data-processing workflow output
    folder (WGX: contains metaphlan/humann/kneaddata output files;
    16S: contains OTU or ASV table).

    Single project:
      nextflow run main.nf --input /path/to/wmgx_output [options]

    Multiple projects (parallel, glob pattern):
      nextflow run main.nf --input "/data/*/wmgx_output" [options]
    ═══════════════════════════════════════════════════════════════
    """.stripIndent()
    System.exit(1)
}

if (!['pdf', 'html'].contains(params.format)) {
    log.error "ERROR: --format must be 'pdf' or 'html'. Got: '${params.format}'"
    System.exit(1)
}

// ============================================================================
// Import modules
// ============================================================================
include { GENERATE_VIS_REPORT } from './modules/generate_report'

// ============================================================================
// Main workflow
// ============================================================================
workflow {

    log.info """
    ┌───────────────────────────────────────────────────────────┐
    │     bioBakery Visualization Workflow  (Nextflow DSL2)      │
    └───────────────────────────────────────────────────────────┘
      Input dir(s)  : ${params.input}
      Output dir    : ${params.output}
      Report format : ${params.format}
      Metadata file : ${params.input_metadata ?: 'none'}
      Profile       : ${workflow.profile ?: 'standard (local)'}
      Container     : ${workflow.containerEngine ?: 'none'}
      Work dir      : ${workflow.workDir}
    ───────────────────────────────────────────────────────────
    """.stripIndent()

    // -------------------------------------------------------------------------
    // Input channel — one item per project directory.
    // Accepts a plain path or a glob pattern, e.g. "/data/*/wmgx_output".
    // -------------------------------------------------------------------------
    Channel
        .fromPath(params.input, type: 'dir', checkIfExists: true)
        .set { input_ch }

    // -------------------------------------------------------------------------
    // Metadata — passed as an absolute-path string (val) so it works on both
    // shared HPC filesystems and cloud storage paths without double-staging.
    // When --input-metadata is not set, the empty string "" is passed and the
    // process omits the flag from the vis.py call.
    // -------------------------------------------------------------------------
    def meta_abs = params.input_metadata
        ? file(params.input_metadata).toAbsolutePath().toString()
        : ""

    if (params.input_metadata && !file(params.input_metadata).exists()) {
        error "ERROR: Metadata file not found: ${params.input_metadata}"
    }

    metadata_val = Channel.value(meta_abs)

    // -------------------------------------------------------------------------
    // Run visualization for each input directory (embarrassingly parallel)
    // -------------------------------------------------------------------------
    GENERATE_VIS_REPORT(
        input_ch,
        metadata_val
    )

    // -------------------------------------------------------------------------
    // Surface output paths in the log
    // -------------------------------------------------------------------------
    GENERATE_VIS_REPORT.out.report
        .ifEmpty { log.warn "No report file(s) found — check process logs." }
        .view    { f -> "  ✓ Report  : ${f}" }

    GENERATE_VIS_REPORT.out.archive
        .view { f -> "  ✓ Archive : ${f}" }
}

// ============================================================================
// Completion / error hooks
// ============================================================================
workflow.onComplete {
    log.info """
    ───────────────────────────────────────────────────────────
    Workflow status : ${workflow.success ? 'COMPLETED ✓' : 'FAILED ✗'}
    Duration        : ${workflow.duration}
    Exit code       : ${workflow.exitStatus}
    Results         : ${params.output}/
    Trace / report  : ${params.output}/pipeline_info/
    ───────────────────────────────────────────────────────────
    """.stripIndent()
}

workflow.onError {
    log.error "Workflow failed: ${workflow.errorMessage}\nSee .nextflow.log for details."
}
