CHANGES IN VERSION 0.99.37
------------------------

NEW FEATURES

    o [09-08-2026] Added GENCODEHub annotation through AnnotationHub and
      rebuilt the Quarto report dashboard with contextual figure and
      workbook notes.

    o [08-08-2026] Added SVA, GLM, LME, Manhattan, Venn/intersection, and
      model-variable figures with selectable report browsers and ENmix controls
      integrated into Quality Control.

    o [08-08-2026] GLM and LME models now support optional coefficient and
      omnibus phenotype Venn outputs with expanded factor terms, positional
      labels, UCSC and GENCODE v50 figures, and ordered workbook sheets.

    o [07-08-2026] Expanded GLM, lmerTest/lme4, and nlme visualisations and the
      Quarto report with model-design and distribution diagnostics,
      annotated-result Manhattan and optional model-level Venn plots,
      selectable workbook sheets, and right-aligned navigation.

    o [05-08-2026] GLM analyses now support optional omnibus F tests for
      phenotype and phenotype-by-interaction terms.

    o [05-08-2026] The Makefile configuration now uses model-specific
      covariate, factor, and scaling variables for GLM and LME analyses.

    o [05-08-2026] Added modelSections to dnamReport() for preprocessing-only,
      GLM-only, LME-only, or complete reports. The complete report remains the
      default.

CHANGES IN VERSION 0.99.36
------------------------

NEW FEATURES

    o [23-07-2026] Reduced GLM, lmerTest/lme4, and nlme EWAS memory use with
      bounded response-block parallelism, memory-aware workers, and resumable
      compact phenotype summaries instead of genome-wide full-model RDS files.

    o [23-07-2026] Recorded conditions in one Model.Message column and reported
      unavailable results through workbook metadata.

    o [22-07-2026] Added lmerTest omnibus F tests for complete
      phenotype main effects or phenotype-by-interaction terms, with
      Satterthwaite or Kenward-Roger denominator degrees of freedom and
      CpG-adjusted results in annotated LME workbooks and reports.

    o [19-07-2026] Added scaleVars predictor standardisation, consistent
      GLM/lme4/nlme metadata, single-file SVA phenotype updates,
      ordered annotated workbooks, and source-defined report card titles.

    o [18-07-2026] Added optional sex-mismatch removal.

    o [14-07-2026] Improved dnamReport() with lightweight paged and filterable
      tables, faster tab loading, and phenotype-labelled GLM, LME, and nlme
      formula notes, with biological participant detection in the Data summary.

    o [09-06-2026] Added methylation-scale support for beta, M-value, and copy
      number phenotype inputs while retaining beta values for cell-composition
      estimation and clock-foundation inputs.

    o [09-06-2026] Improved methylationGLM and methylationLME throughput with
      backend-aware CpG batching, Linux fork support, PSOCK fallback, and
      fit-time summary caching while preserving the glm2 and lmerTest/lme4 model
      engines.

    o [09-06-2026] Added optional nlme-backed methylationLME fitting through
      LME_LIBS/lmeLibs, with none, AR1, and CAR1 residual correlation choices
      while preserving existing LME inputs, output files, and result classes.

    o [09-06-2026] Added LME_CORRELATION_VAR/correlationVar so nlme AR1 and
      CAR1 residual structures can use an explicit within-person ordering
      variable.

    o [09-06-2026] Updated methylationLME fixed-effect assembly so longitudinal
      models are built from phenotypes, covariates, phenotype-specific PRS
      terms, and optional interactions.

    o [07-06-2026] Made the svaEnmix matrix plot adapt to larger
      surrogate-variable matrices, paginate oversized matrices, and suppress
      oversized SentrixID legends with an explicit log note.

    o [07-06-2026] Renamed cross-reactive probe inputs to probe-exclusion inputs,
      with support for multiple files and optional EPICv2 manifest flags.

    o [06-06-2026] Added IDAT_FORCE to optionally force minfi IDAT parsing for
      validated mixed-size IDAT inputs.

    o [05-06-2026] Added configurable cross-reactive probe ID detection. 

    o [30-05-2026] Added an overview vignette with visual summaries of the main
      dnaEPICO functions.

    o [23-05-2026] Added a Quarto dashboard report workflow for dnamReport(). 

    o [14-04-2026] Introduced a modular and reproducible pipeline for preprocessing
      Illumina DNA methylation array data (EPICv2, EPIC and 450K). 

    o [18-02-2026] Started development of the dnaEPICO package.
