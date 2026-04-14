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
#' @param ref Character. Reference panel name. Supported values are `"saliva"`
#'   and `"salivaEPIC"`.
#' @param constrained Logical. If `TRUE`, estimated cell proportions are
#'   constrained to sum to one.
#'
#' @return A `data.table` with one row per sample and one column per estimated
#'   cell type.
#'
#' @references
#' Murat K, et al. Ewastools: Infinium Human Methylation BeadChip pipeline for population epigenetics integrated into Galaxy. *GigaScience*. 2020;9(5):giaa049.
#' Houseman EA, Accomando WP, Koestler DC, et al. DNA methylation arrays as surrogate measures of cell mixture distribution. *BMC Bioinformatics*. 2012;13:86.
#' Reinius LE, Acevedo N, Joerink M, et al. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. *PLoS One*. 2012;7(7):e41361.
#' Bakulski KM, Feinberg JI, Andrews SV, et al. DNA methylation of cord blood cell types: applications for mixed cell birth studies. *Epigenetics*. 2016;11(5):354-362.
#' de Goede OM, Razzaghian HR, Price EM, et al. Nucleated red blood cells impact DNA methylation and expression analyses of cord blood hematopoietic cells. *Clinical Epigenetics*. 2015;7:95.
#' Gervin K, Salas LA, Bakulski KM, et al. Cell type specific DNA methylation in cord blood: a 450K reference data set and cell count-based validation of estimated cell type composition. *Epigenetics*. 2016;11(9):690-698.
#' Gervin K, Salas LA, Bakulski KM, et al. Systematic evaluation and validation of reference and library selection methods for deconvolution of cord blood DNA methylation data. *bioRxiv*. 2019. doi:10.1101/570457.
#' Salas LA, Koestler DC, Butler RA, et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. *Genome Biology*. 2018;19:64.
#' Heiss JA, Just AC, Brenner H. Training a model for estimating leukocyte composition using whole-blood DNA methylation and cell counts as reference. *Epigenomics*. 2017;9(1):13-20.
#' Middleton LYM, Dou J, Mill J, et al. Saliva cell type DNA methylation reference panel for epidemiology studies in children. 2020.

#' @examples
#' ref_file <- system.file("extdata", "saliva.txt", package = "dnaEPICO")
#' ref_panel <- as.matrix(utils::read.table(ref_file))
#' meth <- ref_panel[1:20, , drop = FALSE]
#' colnames(meth) <- c("sample1", "sample2")
#' estimateLC(
#'   meth = meth,
#'   ref = "saliva",
#'   constrained = FALSE
#' )
#'
#' @export
estimateLC <- function(meth, ref, constrained = FALSE) {

    J <- ncol(meth)

    coefs <- utils::read.table(system.file("extdata", paste0(ref, ".txt"), package = "dnaEPICO"))

    coefs <- as.matrix(coefs)
    n_celltypes <- ncol(coefs)


    detect_probe_id_type <- function(x, label) {
      legacy_regex <- "^cg[0-9]{8}$|^rs[0-9]+$|^ch\\.[[:alnum:]_]+\\.\\d+[FR]$"
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

    # Check that row names of coefs and meth are compatible 
    find_matching_rows <- function(query, rownames) {
        if (query_type == row_type) {

          return(match(query, rownames))

        } else if (query_type == "legacy" & row_type == "epicv2") {
          row_loci <- sub("_[TB][CO]\\d+$", "", rownames)

          return(match(query, row_loci))

        } else if (query_type == "epicv2" & row_type == "legacy") {

          stop("Query contains EPICv2 probe IDs but dataset is of legacy type")

        } else {

          stop("This should be a dead branch")
        }
    }

    markers <- find_matching_rows(rownames(coefs), rownames(meth))

    if (any(is.na(markers))) {
      coefs <- coefs[!is.na(markers), , drop = FALSE]
      markers <- stats::na.omit(markers)
    }

    EST <- vapply(seq_len(J), function(j) {

      tmp <- meth[markers, j]
      i <- !is.na(tmp)

      if (constrained == FALSE) {

        quadprog::solve.QP(
          t(coefs[i, ]) %*% coefs[i, ],
          t(coefs[i, ]) %*% tmp[i],
          diag(n_celltypes),
          rep(0, n_celltypes)
        )$solution

      } else {

        quadprog::solve.QP(
          t(coefs[i, ]) %*% coefs[i, ],
          t(coefs[i, ]) %*% tmp[i],
          cbind(1, diag(n_celltypes)),
          c(1, rep(0, n_celltypes)),
          meq = 1
        )$solution
      }
    }, FUN.VALUE = numeric(n_celltypes))

    EST <- t(EST)
    colnames(EST) <- colnames(coefs)
    EST <- data.table::data.table(EST)

    return(EST)
}
