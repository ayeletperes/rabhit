# RAbHIT: R Antibody Haplotype Inference Tool

<!-- badges: start -->
[![R-CMD-check](https://github.com/ayeletperes/rabhit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ayeletperes/rabhit/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/ayeletperes/rabhit/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/ayeletperes/rabhit/actions/workflows/test-coverage.yaml)
[![docs](https://github.com/ayeletperes/rabhit/actions/workflows/docs.yaml/badge.svg)](https://ayeletperes.github.io/rabhit/)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)
<!-- badges: end -->

Analysis of antibody repertoires by high throughput sequencing is of major importance in understanding adaptive immune responses. Our knowledge of variations in the genomic loci encoding antibody genes is incomplete, mostly due to technical difficulties in aligning short reads to these highly repetitive loci. The partial knowledge results in conflicting V-D-J gene assignments between different algorithms, and biased genotype and haplotype inference. Previous studies have shown that haplotypes can be inferred by taking advantage of IGHJ6 heterozygosity, observed in approximately one third of the population.

**RAbHIT is a haplotype infrence tool based on a robust novel method for determining V-D-J haplotypes by adapting a Bayesian framework**. Our method extends haplotype inference to IGHD, IGHV, IGKJ, IGKV, and IGLV based analysis, thereby enabling inference of complex genetic events like deletions and copy number variations in the entire population. Based on this method we developed an R package, which implements the method on sequences from naive B-cells, for both the heavy and the light chains. The package offers a haplotype and single chromosome deletion inference based on an anchor gene.  The inferred haplotypes and deletion patterns may have clinical implications for genetic predispositions to diseases. 


## Core Abilities ##

* Haplotype inference
* Single chromosome deletion detection
* Two chromosome deletion detection

## Required Input ##

* Pre-processed antibody repertoire sequencing data with heterozygosity in at least one gene. Antibody repertoire sequencing data is in a data frame format. Each row represents a unique observation and columns represent data about that observation. The names of the required columns are provided below along with a short description.
* Database of germline gene sequences

| Column name   | Description                                                       |
| ------------- |-------------------------------------------------------------------|
| Subject name  | Subject name                                                      |
| V_CALL        | (Comma separated) name(s) of the nearest V allele(s) (IMGT format)|
| D_CALL        | (Comma separated) name(s) of the nearest D allele(s) (IMGT format)|
| J_CALL        | (Comma separated) name(s) of the nearest J allele(s) (IMGT format)|



## Installation ##

RAbHIT is available for installion either from CRAN or from the development version.

### RAbHIT CRAN installation ###

```R
install.packages("rabhit")
```

### RAbHIT development version installation ###

To install the latest development version directly from GitHub:

```R
# install.packages("devtools")
devtools::install_github("ayeletperes/rabhit")
```

To build from the source code, first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "plotly"))
```

Then clone the repository and build:

```R
library(devtools)
install_deps()
document()
build()
install()
```

## Documentation ##

A complete documentation of RAbHIT is available at: https://ayeletperes.github.io/rabhit/ or in your local repository at: ./vignettes/RAbHIT-vignette.html


## Citation ##

If you use RAbHIT, please cite both of the following papers (you can also run `citation("rabhit")` in R):

> Gidoni M, Snir O, Peres A, *et al.* Mosaic deletion patterns of the human antibody heavy chain gene locus shown by Bayesian haplotyping. *Nature Communications* 10, 628 (2019). doi:[10.1038/s41467-019-08489-3](https://doi.org/10.1038/s41467-019-08489-3)

> Peres A, Gidoni M, Polak P, Yaari G. RAbHIT: R Antibody Haplotype Inference Tool. *Bioinformatics* 35(22):4840–4842 (2019). doi:[10.1093/bioinformatics/btz481](https://doi.org/10.1093/bioinformatics/btz481)

<details>
<summary>BibTeX</summary>

```bibtex
@article{gidoni2019mosaic,
  title     = {Mosaic deletion patterns of the human antibody heavy chain gene locus shown by Bayesian haplotyping},
  author    = {Gidoni, Moriah and Snir, Omri and Peres, Ayelet and Polak, Pazit and Lindeman, Ida and Mikocziova, Ivana and Sarna, Vikas Kumar and Lundin, Knut E. A. and Clouser, Christopher and Vigneault, Francois and Yaari, Gur},
  journal   = {Nature Communications},
  volume    = {10},
  number    = {1},
  pages     = {628},
  year      = {2019},
  doi       = {10.1038/s41467-019-08489-3},
  publisher = {Nature Publishing Group}
}

@article{peres2019rabhit,
  title     = {{RAbHIT}: {R} Antibody Haplotype Inference Tool},
  author    = {Peres, Ayelet and Gidoni, Moriah and Polak, Pazit and Yaari, Gur},
  journal   = {Bioinformatics},
  volume    = {35},
  number    = {22},
  pages     = {4840--4842},
  year      = {2019},
  doi       = {10.1093/bioinformatics/btz481},
  publisher = {Oxford University Press}
}
```

</details>


## Contact ##

For help, questions, or suggestions, please contact:

* [Ayelet Peres](mailto:ayelet.peres@yale.edu)
* [Moriah Gidoni](mailto:moriah.cohen@biu.ac.il)
* [Gur Yaari](mailto:gur.yaari@biu.ac.il)
* [Issue tracker](https://github.com/ayeletperes/rabhit/issues)


## Copying ##

RAbHIT is free for use under the [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode)
