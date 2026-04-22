/*
 * GENERATE_VIS_REPORT
 *
 * Runs the biobakery_workflows vis visualization pipeline.
 * All analysis (alpha diversity, EC annotation, Pweave report generation,
 * and archiving) is performed inside a single process, matching the
 * original AnADAMA2 workflow's task graph exactly.
 *
 * Inputs
 * ------
 * input_dir     : directory of biobakery data-processing workflow output
 *                 (contains metaphlan / humann / kneaddata files for WGX, or
 *                  OTU / ASV tables for 16S)
 * metadata_path : absolute path to optional metadata TSV file, or ""
 *                 (pass as val so it works on shared filesystems + cloud)
 *
 * Outputs
 * -------
 * report     : PDF or HTML report file
 * archive    : ZIP archive containing report + figures + data
 * figures    : figures/ subdirectory (for downstream use / archiving)
 * data_files : data/ subdirectory (TSV tables)
 */

process GENERATE_VIS_REPORT {

    tag { input_dir.name }

    // Publish outputs to: ${params.output}/<project-dir-name>/
    publishDir "${params.output}/${input_dir.name}", mode: 'copy', overwrite: true

    // Retry once on transient failures (LaTeX/Pandoc OOM, network glitch on cloud)
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 143] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path input_dir
    val  metadata_path   // absolute path string or "" (empty = no metadata)

    output:
    path("vis_out/*.${params.format}"),  emit: report,     optional: true
    path("vis_out.zip"),                  emit: archive
    path("vis_out/figures"),              emit: figures,    optional: true
    path("vis_out/data"),                 emit: data_files, optional: true

    script:
    // -------------------------------------------------------------------------
    // Build optional argument strings.
    // Each optional param maps 1-to-1 to the corresponding vis.py flag.
    // List params use repeated flags, matching anadama2's action="append".
    // -------------------------------------------------------------------------
    def meta_arg       = metadata_path
                         ? "--input-metadata '${metadata_path}'"           : ""
    def project_arg    = params.project_name
                         ? "--project-name '${params.project_name}'"       : ""
    def author_arg     = params.author_name
                         ? "--author-name '${params.author_name}'"         : ""
    def header_arg     = params.header_image
                         ? "--header-image '${params.header_image}'"       : ""
    def intro_arg      = params.introduction_text
                         ? "--introduction-text '${params.introduction_text}'" : ""
    def template_arg   = params.use_template
                         ? "--use-template '${params.use_template}'"       : ""
    def picard_arg     = params.input_picard
                         ? "--input-picard '${params.input_picard}'"       : ""
    def picard_ext_arg = "--input-picard-extension '${params.input_picard_extension}'"

    // Expand list parameters into repeated flags
    def cat_args  = (params.metadata_categorical instanceof List
                      ? params.metadata_categorical
                      : params.metadata_categorical ? [params.metadata_categorical] : [])
                     .collect { v -> "--metadata-categorical '${v}'" }.join(" ")

    def cont_args = (params.metadata_continuous instanceof List
                      ? params.metadata_continuous
                      : params.metadata_continuous ? [params.metadata_continuous] : [])
                     .collect { v -> "--metadata-continuous '${v}'" }.join(" ")

    def excl_args = (params.metadata_exclude instanceof List
                      ? params.metadata_exclude
                      : params.metadata_exclude ? [params.metadata_exclude] : [])
                     .collect { v -> "--metadata-exclude '${v}'" }.join(" ")

    def ftype_args = (params.input_file_type instanceof List
                       ? params.input_file_type
                       : params.input_file_type ? [params.input_file_type] : [])
                      .collect { v -> "--input-file-type '${v}'" }.join(" ")

    """
    set -euo pipefail

    # ── Locate the biobakery_workflows vis entry point ─────────────────────
    # 1. The 'biobakery_workflows' CLI (installed via pip / conda)
    # 2. The vis.py script (installed as a script via setup.py)
    # 3. Direct Python module invocation as fallback
    if command -v biobakery_workflows &>/dev/null; then
        VIS_CMD="biobakery_workflows vis"
    elif command -v vis.py &>/dev/null; then
        VIS_CMD="vis.py"
    else
        echo "ERROR: biobakery_workflows not found. Install it or use a container." >&2
        exit 1
    fi

    # ── Absolute path to the staged input directory ────────────────────────
    INPUT_ABS=\$(realpath "${input_dir}")

    # ── Run the visualization workflow ─────────────────────────────────────
    # --local-jobs 1  forces anadama2 to run tasks sequentially within this
    #                 process; prevents nested job scheduler submissions when
    #                 running inside SLURM / PBS / AWS Batch / GCP Batch.
    # --output vis_out keeps output isolated to this Nextflow work directory.
    \${VIS_CMD} \\
        --input  "\${INPUT_ABS}" \\
        --output vis_out \\
        --format ${params.format} \\
        --min-abundance ${params.min_abundance} \\
        --min-samples ${params.min_samples} \\
        --max-sets-heatmap ${params.max_sets_heatmap} \\
        --max-missing ${params.max_missing} \\
        --max-sets-barplot ${params.max_sets_barplot} \\
        --max-groups-barplot ${params.max_groups_barplot} \\
        --correlation-threshold ${params.correlation_threshold} \\
        ${meta_arg} \\
        ${project_arg} \\
        ${author_arg} \\
        ${header_arg} \\
        ${intro_arg} \\
        ${template_arg} \\
        ${picard_arg} \\
        ${picard_ext_arg} \\
        ${cat_args} \\
        ${cont_args} \\
        ${excl_args} \\
        ${ftype_args} \\
        --local-jobs 1

    # ── Verify expected outputs are present ────────────────────────────────
    if [ ! -f vis_out.zip ]; then
        echo "ERROR: vis_out.zip was not created. Check logs above." >&2
        exit 1
    fi

    # Summarise what was produced
    REPORT_COUNT=\$(find vis_out/ -maxdepth 1 -name "*.${params.format}" 2>/dev/null | wc -l || echo 0)
    FIGURE_COUNT=\$(find vis_out/figures/ -name "*.png" 2>/dev/null | wc -l || echo 0)
    echo "────────────────────────────────────────────────────────"
    echo "Project  : ${input_dir.name}"
    echo "Format   : ${params.format}"
    echo "Reports  : \${REPORT_COUNT}"
    echo "Figures  : \${FIGURE_COUNT} PNG files"
    echo "Archive  : vis_out.zip"
    echo "────────────────────────────────────────────────────────"
    """
}
