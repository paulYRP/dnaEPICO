#' Estimate saliva cell proportions from DNA methylation beta values
#'
#' Estimate cell-type proportions with the saliva reference panels bundled in
#' `dnaEPICO`. This function keeps the original `estimateLC()` interface used by
#' the package while using the internal reference files distributed in
#' `inst/extdata`.
#'
#' @param meth Numeric matrix of beta values with CpGs in rows and samples in
#'   columns. Row names must contain probe identifiers compatible with the
#'   selected reference.
#' @param ref Character. Reference panel name. Supported values are `'saliva'`
#'   and `'salivaEPIC'`.
#' @param constrained Logical. If `TRUE`, estimated cell proportions are
#'   constrained to sum to one.
#'
#' @return A `data.table` with one row per sample and one column per estimated
#'   cell type.
#'
#' @references
#' Murat K, et al. Ewastools: Infinium Human Methylation BeadChip pipeline for
#' population epigenetics integrated into Galaxy. *GigaScience*.
#' 2020;9(5):giaa049.
#' Houseman EA, Accomando WP, Koestler DC, et al. DNA methylation arrays as
#' surrogate measures of cell mixture distribution. *BMC Bioinformatics*.
#' 2012;13:86.
#' Reinius LE, Acevedo N, Joerink M, et al. Differential DNA methylation in
#' purified human blood cells: implications for cell lineage and studies on
#' disease susceptibility. *PLoS One*. 2012;7(7):e41361.
#' Bakulski KM, Feinberg JI, Andrews SV, et al. DNA methylation of cord blood
#' cell types: applications for mixed cell birth studies. *Epigenetics*.
#' 2016;11(5):354-362.
#' de Goede OM, Razzaghian HR, Price EM, et al. Nucleated red blood cells
#' impact DNA methylation and expression analyses of cord blood hematopoietic
#' cells. *Clinical Epigenetics*. 2015;7:95.
#' Gervin K, Salas LA, Bakulski KM, et al. Cell type specific DNA methylation
#' in cord blood: a 450K reference data set and cell count-based validation of
#' estimated cell type composition. *Epigenetics*. 2016;11(9):690-698.
#' Gervin K, Salas LA, Bakulski KM, et al. Systematic evaluation and validation
#' of reference and library selection methods for deconvolution of cord blood
#' DNA methylation data. *bioRxiv*. 2019. doi:10.1101/570457.
#' Salas LA, Koestler DC, Butler RA, et al. An optimized library for
#' reference-based deconvolution of whole-blood biospecimens assayed using the
#' Illumina HumanMethylationEPIC BeadArray. *Genome Biology*. 2018;19:64.
#' Heiss JA, Just AC, Brenner H. Training a model for estimating leukocyte
#' composition using whole-blood DNA methylation and cell counts as reference.
#' *Epigenomics*. 2017;9(1):13-20.
#' Middleton LYM, Dou J, Mill J, et al. Saliva cell type DNA methylation
#' reference panel for epidemiology studies in children. 2020.

#' @examples
#' ref_file <- system.file("extdata", "saliva.txt", package = "dnaEPICO")
#' ref_panel <- as.matrix(utils::read.table(ref_file))
#' meth <- ref_panel[1:20, , drop = FALSE]
#' colnames(meth) <- c("sample1", "sample2")
#' estimateLC(
#'     meth = meth,
#'     ref = "saliva",
#'     constrained = FALSE
#' )
#'
#' @export
estimateLC <- function(meth, ref, constrained = FALSE) {
    constrained <- validateLogicalScalarDnaEpico(
        constrained,
        "constrained"
    )
    validateCellCompositionBetaStructureDnaEpico(
        beta = meth,
        requireSampleNames = FALSE, objectName = "meth"
    )
    validateCellCompositionReferenceDnaEpico(ref)
    methylation_range <- summarizeMethylationRangeDnaEpico(meth, "beta")
    estimate <- estimateLCFromBetaDnaEpico(
        meth = meth,
        ref = ref, constrained = constrained
    )
    attr(estimate, "methylationRange") <- methylation_range
    estimate
}

validateCellCompositionBetaStructureDnaEpico <- function(
    beta,
    requireSampleNames = TRUE, objectName = "beta"
) {
    if (!is.matrix(beta) || !is.numeric(beta) || nrow(beta) ==
        0L || ncol(beta) == 0L) {
        stop(objectName, " must be a non-empty numeric matrix.",
            call. = FALSE
        )
    }
    probe_ids <- rownames(beta)
    if (is.null(probe_ids)) {
        stop(objectName, " must have probe identifiers in row names.",
            call. = FALSE
        )
    }
    validateMethylationProbeIdentifiersDnaEpico(probe_ids, paste(
        objectName,
        "row names"
    ))
    if (isTRUE(requireSampleNames) && is.null(colnames(beta))) {
        stop(objectName, " must have sample names in columns.",
            call. = FALSE
        )
    }
    if (!is.null(colnames(beta))) {
        validateSampleIdentifiersDnaEpico(colnames(beta), paste(
            objectName,
            "sample identifiers"
        ))
    }
    invisible(TRUE)
}

validateCellCompositionReferenceDnaEpico <- function(ref) {
    if (!is.character(ref) || length(ref) != 1L || is.na(ref) ||
        !(ref %in% c("saliva", "salivaEPIC"))) {
        stop("Unsupported ref. Supported values are 'saliva' and 'salivaEPIC'.",
            call. = FALSE
        )
    }
    invisible(ref)
}

estimateLCFromBetaDnaEpico <- function(meth, ref, constrained) {
    ref_file <- system.file("extdata", paste0(ref, ".txt"),
        package = "dnaEPICO")
    if (!nzchar(ref_file)) {
        stop("The bundled reference file was not found for ref '",
            ref, "'.",
            call. = FALSE
        )
    }

    J <- ncol(meth)

    coefs <- utils::read.table(ref_file)

    coefs <- as.matrix(coefs)
    n_celltypes <- ncol(coefs)


    detect_probe_id_type <- function(x, label) {
        legacy_regex <-
            "^cg[0-9]{8}$|^rs[0-9]+$|^ch\\.[[:alnum:]_]+\\.\\d+[FR]$"
        epicv2_regex <- "^(cg|ch|rs|nv).+_[TB][CO][0-9]+$"

        if (all(grepl(legacy_regex, x))) {
            return("legacy")
        }

        if (all(grepl(epicv2_regex, x))) {
            return("epicv2")
        }

        stop(label, " is (partly) invalid!")
    }

    query_type <- detect_probe_id_type(rownames(coefs), "Query")
    row_type <- detect_probe_id_type(rownames(meth), "Input row names")

    if (identical(query_type, row_type)) {
        markers <- match(rownames(coefs), rownames(meth))
        retained_markers <- !is.na(markers)
        coefs <- coefs[retained_markers, , drop = FALSE]
        marker_meth <- meth[markers[retained_markers], , drop = FALSE]
    } else if (query_type == "legacy" && row_type == "epicv2") {
        row_loci <- sub("_[TB][CO]\\d+$", "", rownames(meth))
        marker_groups <- lapply(rownames(coefs), function(marker) {
            which(row_loci ==
                marker)
        })
        retained_markers <- lengths(marker_groups) > 0L
        coefs <- coefs[retained_markers, , drop = FALSE]
        marker_groups <- marker_groups[retained_markers]
        marker_meth <- t(vapply(marker_groups, function(indices) {
            locus_values <- meth[indices, , drop = FALSE]
            observed <- colSums(!is.na(locus_values))
            locus_means <- rep(NA_real_, J)
            has_values <- observed > 0L
            locus_means[has_values] <- colSums(locus_values[,
                has_values,
                drop = FALSE
            ], na.rm = TRUE) / observed[has_values]
            locus_means
        }, numeric(J)))
    } else if (query_type == "epicv2" && row_type == "legacy") {
        stop("Query contains EPICv2 probe IDs but dataset is of legacy type")
    } else {
        stop("Unsupported probe identifier combination.", call. = FALSE)
    }

    if (nrow(marker_meth) < n_celltypes || qr(coefs)$rank < n_celltypes) {
        stop("Too few independent reference markers overlap meth for cell ",
            "composition estimation.",
            call. = FALSE
        )
    }

    EST <- vapply(seq_len(J), function(j) {
        tmp <- marker_meth[, j]
        i <- !is.na(tmp)

        if (sum(i) < n_celltypes || qr(coefs[i, , drop = FALSE])$rank <
            n_celltypes) {
            stop("Sample ", if (is.null(colnames(meth))) {
                j
            } else {
                colnames(meth)[[j]]
            }, " has too few independent markers for cell composition estimation.",
            call. = FALSE
            )
        }

        if (!isTRUE(constrained)) {
            quadprog::solve.QP(
                t(coefs[i, ]) %*% coefs[i, ],
                t(coefs[i, ]) %*% tmp[i], diag(n_celltypes),
                rep(0, n_celltypes)
            )$solution
        } else {
            quadprog::solve.QP(t(coefs[i, ]) %*% coefs[i, ],
                t(coefs[i, ]) %*% tmp[i], cbind(1, diag(n_celltypes)),
                c(1, rep(0, n_celltypes)),
                meq = 1
            )$solution
        }
    }, FUN.VALUE = numeric(n_celltypes))

    EST <- t(EST)
    colnames(EST) <- colnames(coefs)
    EST <- data.table::data.table(EST)

    EST
}
