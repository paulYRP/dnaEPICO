#' Run estimateLC.R
#' @import stringr
#' @importFrom data.table data.table setkey
#' @import quadprog
#'
#' @param meth Matrix of beta values
#' @param ref Choice of reference dataset: available options are `saliva` and `salivaEPIC`. The functions are extracted from `ewastools (Murat, K, et al. 2020)` and adapted to work onloy with the saliva reference panel of `Middleton et al. 2020`.
#' @param constrained Force that all cell proportions sum up to 1.
#'
#' @return A data.table containing estimated cell proportions for each sample. 
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
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' estimateLC(
#'  meth = beta, 
#'  ref = "saliva",
#'  constrained = FALSE
#' )
#' }
#'
#' @export
estimateLC <- function(meth, ref, constrained = FALSE) {

  J <- ncol(meth)

  coefs <- read.table(system.file("extdata", paste0(ref, ".txt"), package = "dnaEPICO"))

  coefs <- as.matrix(coefs)
  n_celltypes <- ncol(coefs)

  # Check that row names of coefs and meth are compatible 
  find_matching_rows <- function(query, rownames) {
      legacy_regex <- stringr::regex("^cg\\d{8}$|^rs\\d+$|^ch\\.\\w+\\.\\d+[FR]$")
      epicv2_regex <- stringr::regex("^[cg|ch|rs|nv].*_[TB][CO]\\d+$")

      if (all(stringr::str_detect(query, pattern = legacy_regex))) {
        query_type <- "legacy"
      } else if (all(stringr::str_detect(query, pattern = epicv2_regex))) {
        query_type <- "epicv2"
      } else {
        stop("Query is (partly) invalid!")
      }

      if (all(stringr::str_detect(rownames, pattern = legacy_regex))) {
        row_type <- "legacy"
      } else if (all(stringr::str_detect(rownames, pattern = epicv2_regex))) {
        row_type <- "epicv2"
      } else {
        stop("(Some) row names are invalid!")
      }

      if (query_type == row_type) {

        return(match(query, rownames))

      } else if (query_type == "legacy" & row_type == "epicv2") {

        df_map <- data.table::data.table(
          MANIFESTS$EPICv2[, c("ilmn_id", "probe_id")]
        )

        data.table::setkey(df_map, probe_id)

        query <- df_map[query, nomatch = NA, mult = "first"]$ilmn_id

        return(match(query, rownames))

      } else if (query_type == "epicv2" & row_type == "legacy") {

        stop("Query contains EPICv2 probe IDs but dataset is of legacy type")

      } else {

        stop("This should be a dead branch")
      }
  }

  markers <- find_matching_rows(rownames(coefs), rownames(meth))

  if (any(is.na(markers))) {
    coefs <- coefs[!is.na(markers), , drop = FALSE]
    markers <- na.omit(markers)
  }

  EST <- sapply(1:J, function(j) {

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
  })

  EST <- t(EST)
  colnames(EST) <- colnames(coefs)
  EST <- data.table::data.table(EST)

  return(EST)
}
