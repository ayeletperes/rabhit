# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RAbHIT (R Antibody Haplotype Inference Tool) is an R package implementing a Bayesian framework for inferring V-D-J haplotypes and gene deletions from AIRR-seq data. It supports IGH, IGK, IGL, and TRB chains. Published in Gidoni et al. (2019) Nature Communications and Peres & Gidoni et al. (2019) Bioinformatics.

## Build & Development Commands

```bash
# Generate documentation (NAMESPACE, man/ pages) via Roxygen2
Rscript -e "devtools::document()"

# Install package locally
R CMD INSTALL --no-multiarch --with-keep.source --resave-data .

# Run R CMD check (as-cran)
R CMD check --as-cran --timings .

# Build package tarball
R CMD build --resave-data=best --compact-vignettes=both .

# Build vignette
Rscript -e "rmarkdown::render('vignettes/RAbHIT-vignette.Rmd')"

# Regenerate system data (sysdata.rda) from raw files
Rscript data-raw/GenerateSysData.R
```

There is no formal test suite (`tests/` directory). The vignette serves as the integration test.

## Architecture

### Source Files (R/)

- **functions.R** — 7 core exported analysis functions: `createFullHaplotype()` (anchor-based haplotype inference), `geneUsage()`, `deletionsByBinom()` (double chromosome deletions via binomial test), `deletionsByVpooled()` (single chromosome deletions), `nonReliableVGenes()`, `createHaplotypeTable()`, `readHaplotypeDb()`
- **internal_functions.R** — Internal helpers. Key function: `get_probabilites_with_priors()` implements the core Dirichlet-multinomial Bayesian model with epsilon adjustment
- **graphic_functions.R** — 6 exported visualization functions: `plotHaplotype()`, `hapHeatmap()`, `hapDendo()`, `deletionHeatmap()`, `plotDeletionsByBinom()`, `plotDeletionsByVpooled()`. Uses ggplot2 and plotly for interactive HTML5 output
- **Data.R** — Roxygen documentation for bundled datasets (germline references, example data)
- **zzz.R** — Package load message

### Data

- **data/** — Pre-compiled .rda files: germline databases (HVGERM, HDGERM, HJGERM, KVGERM, KJGERM, LVGERM, LJGERM), gene chromosome ordering (GENE.loc), example repertoire (samples_db), example results (samplesHaplotype)
- **data-raw/** — Source FASTA files and binomial test cutoff tables used to generate `R/sysdata.rda` (contains PSEUDO gene lists and Binom.test.gene.cutoff)

### Key Algorithmic Details

- Bayesian haplotype inference uses Dirichlet-multinomial priors comparing H1 (single chromosome) vs H2 (both chromosomes), reported as Bayes factor lK = log10(likelihood ratio)
- Chain type (IGH/IGK/IGL/TRB) determines which germline databases, gene orderings, and binomial cutoffs are used
- `createFullHaplotype()` is the main entry point; it calls `createHaplotypeTable()` internally for each anchor gene

### Dependencies

Core: dplyr, data.table, ggplot2, plotly, alakazam, tigger. Full list in DESCRIPTION Imports. No compiled code (pure R).
