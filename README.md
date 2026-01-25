
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dnaEPICO

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/paulYRP/dnaEPICO)](https://github.com/paulYRP/dnaEPICO/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/paulYRP/dnaEPICO)](https://github.com/paulYRP/dnaEPICO/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- [![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/dnaEPICO.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/dnaEPICO) -->
<!-- [![Bioc devel status](http://www.bioconductor.org/shields/build/devel/bioc/dnaEPICO.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/dnaEPICO) -->
<!-- [![Bioc downloads rank](https://bioconductor.org/shields/downloads/release/dnaEPICO.svg)](http://bioconductor.org/packages/stats/bioc/dnaEPICO/) -->
<!-- [![Bioc support](https://bioconductor.org/shields/posts/dnaEPICO.svg)](https://support.bioconductor.org/tag/dnaEPICO) -->
<!-- [![Bioc history](https://bioconductor.org/shields/years-in-bioc/dnaEPICO.svg)](https://bioconductor.org/packages/release/bioc/html/dnaEPICO.html#since) -->
<!-- [![Bioc last commit](https://bioconductor.org/shields/lastcommit/devel/bioc/dnaEPICO.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/dnaEPICO/) -->
<!-- [![Bioc dependencies](https://bioconductor.org/shields/dependencies/release/dnaEPICO.svg)](https://bioconductor.org/packages/release/bioc/html/dnaEPICO.html#since) -->
[![Codecov test
coverage](https://codecov.io/gh/paulYRP/dnapipeR/graph/badge.svg)](https://app.codecov.io/gh/paulYRP/dnapipeR)
<!-- badges: end -->

The goal of **`dnaEPICO`** is to provide a **modular, reproducible, and
pipeline** for the preprocessing and statistical analysis of Illumina
DNA methylation array data (EPICv2, EPIC and 450K).

The package integrates preprocessing, quality control, phenotype
merging, generalised linear models (GLM), linear mixed-effects models
(LME), and automated report generation. It is designed to run seamlessly
on local machines as well as High-Performance Computing (HPC)
environments via a **GNU Make–based workflow**.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `dnaEPICO` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("dnaEPICO")
```

And the development version from
[GitHub](https://github.com/paulYRP/dnaEPICO) with:

``` r
BiocManager::install("paulYRP/dnaEPICO")
```

## Articles:

- [**A Pilot Epigenome-Wide Study of Posttraumatic Growth: Identifying
  Novel Candidates for Future
  Research**](https://www.mdpi.com/2075-4655/9/4/39)

## Tutorials:

- [**DNA Methylation
  Tutorial**](https://paulYRP.github.io/2025-cpgpneurogenomics-workshop/tutorial.html)
- [**Getting
  Started**](https://github.com/paulYRP/dnaEPICO/wiki/Getting-Started)
- [**Requirements**](https://github.com/paulYRP/dnaEPICO/wiki/Requirements)

## Citation

Below is the citation output from using `citation('dnaEPICO')` in R.
Please run this yourself to check for any updates on how to cite
**dnaEPICO**.

``` r
print(citation('dnaEPICO'), bibtex = TRUE)
#> dnaEPICO: Analysis Pipeline for Illumina EPICv2 DNA Methylation Array
#> Data
#> 
#>   Ruiz P., Metha D. dnaEPICO: Analysis Pipeline for Illumina EPICv2 DNA
#>   Methylation Array Data. Bioconductor package.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {dnaEPICO: Analysis Pipeline for Illumina EPICv2 DNA Methylation Array Data},
#>     author = {Paul Ruiz and Divya Metha},
#>     year = {2026},
#>     url = {https://bioconductor.org/packages/dnaEPICO},
#>   }
```

Please note that the `dnaEPICO` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `dnaEPICO` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.20/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://paulYRP.github.io/dnaEPICO) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.20/biocthis)*.
