CHANGES IN VERSION 0.99.32
------------------------

NEW FEATURES

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

    o [30-05-2026] Added a overview vignette with visual summaries of the main
      dnaEPICO functions.

    o [23-05-2026] Added a Quarto dashboard report workflow for dnamReport(). 

    o [14-04-2026] Introduced a modular and reproducible pipeline for preprocessing
      Illumina DNA methylation array data (EPICv2, EPIC and 450K). 

    o [18-02-2026] Initial Bioconductor submission of the dnaEPICO package. 
