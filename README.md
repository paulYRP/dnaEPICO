
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dnaEPICO <img src="https://raw.githubusercontent.com/paulYRP/dnaEPICO/devel/inst/extdata/dnaEPICORM.svg" align="right" width="140" style="width:140px; height:auto; float:right; margin:0 -24px 16px 20px;" />

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/paulYRP/dnaEPICO)](https://github.com/paulYRP/dnaEPICO/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/paulYRP/dnaEPICO)](https://github.com/paulYRP/dnaEPICO/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![codecov](https://codecov.io/gh/paulYRP/dnaEPICO/graph/badge.svg?token=5ZD6K3SMHB)](https://codecov.io/gh/paulYRP/dnaEPICO)

<!-- badges: end -->

<br clear="right" /> <br /> <br />

<p align="center">

<img src="https://raw.githubusercontent.com/paulYRP/dnaEPICO/devel/inst/extdata/dnaEPICO.gif" alt="dnaEPICO workflow preview" width="650" />
</p>

dnaEPICO provides a **modular and reproducible pipeline** for
preprocessing and statistically analysing Illumina DNA methylation array
data from the EPICv2, EPIC, and 450K platforms.

The package supports CpG-wise generalised linear models and longitudinal
mixed-effects models using lmerTest/lme4 or nlme. It also integrates
preprocessing, quality control, phenotype preparation, and automated
reporting for local and high-performance computing (HPC) environments
through a **GNU Make-based workflow**.

Optional omnibus F tests jointly evaluate complete phenotype main
effects or phenotype-by-interaction terms in GLM and lmerTest/lme4
analyses.

## Installation

Install a current `R` release from [CRAN](https://cran.r-project.org/),
then install the development version of `dnaEPICO` from GitHub:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("paulYRP/dnaEPICO")
```

## Articles

- [**A Pilot Epigenome-Wide Study of Posttraumatic Growth: Identifying
  Novel Candidates for Future
  Research**](https://www.mdpi.com/2075-4655/9/4/39)

## Tutorials

- [**DNA Methylation
  Tutorial**](https://paulYRP.github.io/2025-cpgpneurogenomics-workshop/tutorial.html)
- [**Getting
  Started**](https://github.com/paulYRP/dnaEPICO/wiki/Getting-Started)
- [**Requirements**](https://github.com/paulYRP/dnaEPICO/wiki/Requirements)

## Citation

Run `citation("dnaEPICO")` to obtain the current citation:

``` r
print(citation("dnaEPICO"), bibtex = TRUE)
```

    ## To cite dnaEPICO, use:
    ## 
    ##   Ruiz P (2026). "dnaEPICO: Analysis Pipeline for Illumina DNA
    ##   Methylation Array Data." _Epigenomes_. doi:10.3390/epigenomes9040039
    ##   <https://doi.org/10.3390/epigenomes9040039>,
    ##   <https://github.com/paulYRP/dnaEPICO>.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {dnaEPICO: Analysis Pipeline for Illumina DNA Methylation Array Data},
    ##     doi = {10.3390/epigenomes9040039},
    ##     journal = {Epigenomes},
    ##     author = {Paul Ruiz},
    ##     year = {2026},
    ##     url = {https://github.com/paulYRP/dnaEPICO},
    ##   }

`dnaEPICO` builds on R and bioinformatics software cited in the
vignettes and package publications.

## Code of Conduct

The `dnaEPICO` project follows the [Bioconductor Code of
Conduct](https://bioconductor.org/about/code-of-conduct/). Contributors
agree to follow its terms.

## Development tools

- [GitHub Actions](https://github.com/paulYRP/dnaEPICO/actions) runs
  package checks with Bioconductor containers and
  *[BiocCheck](https://bioconductor.org/packages/3.20/BiocCheck)*.
- [Codecov](https://codecov.io/gh/paulYRP/dnaEPICO) and
  *[covr](https://CRAN.R-project.org/package=covr)* report code
  coverage.
- *[pkgdown](https://CRAN.R-project.org/package=pkgdown)* builds the
  [documentation website](https://paulYRP.github.io/dnaEPICO).
- *[styler](https://CRAN.R-project.org/package=styler)* formats R code.
- *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)* generate
  package documentation.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.20/biocthis)*.
