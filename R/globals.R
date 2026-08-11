# ================================================================
# dnaEPICO - Global variable declarations for NSE compatibility
# ================================================================
#
# Purpose:
# Declare variables used via Non-Standard Evaluation (NSE) so that
# R CMD check does not produce "no visible binding for global variable"
#
# This file does not modify runtime behaviour. It only informs the
# static code analyser about variables used inside:
#   - data.table operations
#   - modelling outputs
#   - annotation tables
#
# ================================================================

utils::globalVariables(c(
    # ---------- Illumina annotation / manifest ----------
    "IlmnID", "Name", "AddressA_ID", "AddressB_ID", "Next_Base",
    "Color_Channel", "Probe_Type", "Infinium_Design_Type",
    "CHR", "MAPINFO", "Rep_Num",
    "CHR_37", "MAPINFO_37",
    "CHR_hg38", "Start_hg38", "Strand_hg38",

    # ---------- Chromosome / probe mapping ----------
    "chr37", "chr38", "probe_id",
    "probe_design", "channel", "probe_type", "next_base",

    # ---------- DNAm metrics ----------
    "m", "cn",

    # ---------- Cell composition / estimateLC ----------
    "CpG", "FDR", "P.value", "Std.Error",

    # ---------- GLM / LME outputs ----------
    "Pr(>|t|)", "ResidualSD",
    "value", "expected", "observed", "meanMethylation", "diagnosticP",
    "diagnosticY",

    # ---------- Pipeline variables ----------
    "idatFolder", "rBaseDir"
))

# The report builder runs short stages in a shared state environment. These
# names are supplied by that environment, while the plotting names are ggplot2
# aesthetics. Declaring them here informs static checking only.
utils::globalVariables(base::strsplit(base::paste(
        "additional_downloads analysis annotation",
        "annotation_description annotation_title assets_dir",
        "assets_figures_dir assets_logs_dir available_stage_counts",
        "base_log_cards batch_effect_figure_descriptions",
        "batch_effect_figure_titles batch_effect_notes",
        "batch_effect_page batch_effect_paragraphs batch_items batch_k",
        "browser_id browser_notes build_data_frame_table_section",
        "build_figure_browser build_figure_sections",
        "build_model_formula_notes build_model_visualisation_tabs",
        "build_report_page build_result_table_section",
        "build_workbook_table_section callout_lines chunk_size",
        "collapse_values combine_venn_interaction_tokens compose_page",
        "content_description_html content_note_block copy_figure_assets",
        "copy_first_existing_log_asset copy_log_asset correlation",
        "cpgDetectionPath cpg_detection_path cpg_detection_summary",
        "cpg_summary cpg_table_notes csv_name data data_notes data_page",
        "data_paragraphs data_path data_summary data_table_assets",
        "describe_available_annotations descriptions detPPath",
        "detPThreshold detection_table_warning detection_tables",
        "detp_path dictionary empty_message enmixTab",
        "enmix_figure_descriptions enmix_figure_titles enmix_items",
        "enmix_notes error_message explicit_phenotypes",
        "extract_report_count extract_report_formula_phenotype",
        "extract_sva_log_summary fig_dir figure_data",
        "figure_descriptions figure_items figure_titles",
        "figure_viewer_js_source_path filename find_quarto format_count",
        "format_decimal format_model_formula formulas front_columns",
        "glmTab glm_dir glm_items glm_notes glm_page glm_paragraphs",
        "glm_table_assets glm_table_description glm_table_path",
        "glm_table_title glm_venn_items glm_visualisation_bullets",
        "glm_visualisation_notes glm_visualisations_page",
        "has_lme_interaction html_bullet_list html_escape",
        "html_paragraph html_section id_column imagePattern",
        "image_pattern include_glm include_lme infer_model_from_path",
        "is_blank is_numeric_like items js_quote_result_values",
        "js_result_array known_phenotypes last_log_field lmeTab",
        "lme_analysis_label lme_backend lme_dir lme_interaction_term",
        "lme_items lme_notes lme_page lme_paragraphs lme_table_assets",
        "lme_table_description lme_table_path lme_table_title",
        "lme_venn_items lme_visualisation_bullets",
        "lme_visualisation_notes lme_visualisations_page logTab",
        "log_assets log_capture log_field log_path log_render_helper",
        "log_status_note log_summary logoPath logo_source_path logs",
        "logs_bullets logs_dir logs_notes logs_page magick_available",
        "make_batch_effect_notes make_data_notes make_logs_notes",
        "make_metrics_notes make_quality_control_notes",
        "make_report_observations make_report_overlap_notes",
        "make_report_preprocessing_notes make_report_summary_items",
        "metadata metadata_frame_to_list methylation_stage metricTab",
        "metrics_figure_descriptions metrics_figure_titles",
        "metrics_items metrics_metadata metrics_notes metrics_page",
        "metrics_paragraphs modelSections model_figure_description",
        "model_figure_title model_log_path model_name model_navbar",
        "model_notes model_sections normalization normalize_log_field",
        "observation_bullets outputDir overlap_notes overview",
        "parse_log_fields participant_col person phenoTab pheno_file",
        "pick_column pick_participant_column pick_timepoint_column",
        "plural postprocessing_dir prepare_csv_table_assets",
        "prepare_report_text prepare_table_viewer_directory",
        "prepare_xlsx_table_assets prepare_xlsx_workbook_assets",
        "prepared_report preprocessing_dir preprocessing_notes",
        "preprocessing_paragraphs pretty_label pretty_model_term",
        "primary_sheet projectName project_dir project_name qcTab",
        "qc_dir qc_items quality_control_figure_descriptions",
        "quality_control_figure_titles quality_control_metadata",
        "quality_control_notes quality_control_page",
        "quality_control_paragraphs quality_figure_notes quarto_bin",
        "quarto_yml r_string read_detection_tables",
        "read_optional_workbook_sheet recursive",
        "remove_unrequested_rendered_assets render_status rendered_file",
        "report_inline_markup report_interactions_js_source_path",
        "report_note_text report_number_markup report_observations",
        "report_overlap_notes report_page report_preprocessing_notes",
        "report_sections report_summary_items representative_log_count",
        "resolve_lme_report_metadata resolve_report_formula_records",
        "resolve_report_path result_columns result_dir root_dir",
        "safe_read_table_file same_report_path sampleDetectionPath",
        "sample_detection_path sample_detection_summary sample_summary",
        "sample_table_notes search_index_markup selector_id",
        "sentence_case sentence_list sheet sheet_description",
        "sheet_names sidecar significant site_css slash slugify",
        "sort_values source_data_path source_files source_mode",
        "source_path split_report_metadata_values stage_sample_count",
        "stem subset_figure_items summarize_cpg_detection",
        "summarize_dataset summarize_logs summarize_sample_detection",
        "summarize_sva summary_items survey_columns svaTab sva_dir",
        "sva_summary table_data table_viewer_asset_result terms",
        "text_name threshold_description threshold_title title",
        "title_prefix titles unrequested_table_assets var_prefix",
        "variable venn_effect_label venn_figure_metadata verbose",
        "viewer_assets viewer_js_source_path workbook_assets",
        "workbook_available workbook_href workbook_metadata",
        "write_delimited_table_viewer_assets write_table_viewer_assets",
        "write_table_viewer_chunk write_table_viewer_manifest",
        "write_utf8"
    ), " ", fixed = TRUE)[[1L]])
