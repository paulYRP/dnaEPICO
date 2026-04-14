
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dnaEPICO

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/paulYRP/dnaEPICO)](https://github.com/paulYRP/dnaEPICO/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/paulYRP/dnaEPICO)](https://github.com/paulYRP/dnaEPICO/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![codecov](https://codecov.io/gh/paulYRP/dnaEPICO/graph/badge.svg?token=5ZD6K3SMHB)](https://codecov.io/gh/paulYRP/dnaEPICO)

<!-- badges: end -->

dnaEPICO aims to provide a **modular and reproducible pipeline** for
preprocessing and statistically analysing Illumina DNA methylation array
data from the EPICv2, EPIC, and 450K platforms.

The package integrates preprocessing, quality control, phenotype
merging, generalised linear models (GLM), linear mixed-effects models
(LME), and automated report generation. It is designed to operate on
both local machines and High-Performance Computing (HPC) environments
through a **GNU Make–based workflow**.

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

Alternatively, install directly from GitHub:

``` r
install.packages('devtools')
devtools::install_github('paulYRP/dnaEPICO')
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
print(citation("dnaEPICO"), bibtex = TRUE)
#> dnaEPICO: Analysis Pipeline for Illumina DNA Methylation Array Data,
#> generated as part of 10.3390/epigenomes9040039
#> 
#>   Ruiz P, Mehta D (2025). "dnaEPICO: Analysis Pipeline for Illumina DNA
#>   Methylation Array Data." _Epigenomes_. doi:10.3390/epigenomes9040039
#>   <https://doi.org/10.3390/epigenomes9040039>,
#>   <https://bioconductor.org/packages/dnaEPICO>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {dnaEPICO: Analysis Pipeline for Illumina DNA Methylation Array Data},
#>     doi = {10.3390/epigenomes9040039},
#>     journal = {Epigenomes},
#>     author = {Paul Ruiz and Divya Mehta},
#>     year = {2025},
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
