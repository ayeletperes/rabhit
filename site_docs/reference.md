# Function reference

The complete, always-current documentation for every function and argument is available from R via
`?functionName` and `help(package = "rabhit")`. This page summarizes the exported API.

## Haplotype inference

### `createFullHaplotype()`

Anchor-based V-D-J haplotype inference — the main entry point.

```r
createFullHaplotype(clip_db, toHap_col = c("v_call", "d_call"),
                    hapBy_col = "j_call", hapBy = "IGHJ6", toHap_GERM = NULL,
                    relative_freq_priors = TRUE, kThreshDel = 3, rmPseudo = TRUE,
                    deleted_genes = c(), nonReliable_Vgenes = c(),
                    min_minor_fraction = 0.3, single_gene = TRUE,
                    chain = c("IGH", "IGK", "IGL", "TRB"))
```

Key arguments: `toHap_col` the segment columns to haplotype; `hapBy_col`/`hapBy` the anchor column and
gene; `toHap_GERM` the germline (defaults to the bundled set for `chain`); `kThreshDel` the Bayes-factor
threshold for calling deletions; `min_minor_fraction` the minimum minor-allele fraction for a usable
anchor. Returns a haplotype table with, per gene, the allele assigned to each anchor chromosome, the
candidate alleles, read counts, and the Bayes factor `lK`.

### `createHaplotypeTable()`

Lower-level routine that haplotypes a single gene against the anchor counts using the
Dirichlet-multinomial Bayesian model. Called internally by `createFullHaplotype()` per gene.

### `readHaplotypeDb()`

Read a previously saved haplotype table back into R.

### `convertToASC()`

Prepare data and germline for haplotype inference using PIgLET Allele Similarity
Clusters (ASC), which collapse near-identical/duplicated IGHV alleles and align with
OGRDB reference sets. Returns ASC-named `clip_db`, ASC `germline`, a `genes_order`,
and the `allele_cluster_table`. Requires the optional `piglet` package. See
[Allele clusters (ASC/OGRDB)](allele_clusters.md).

## Gene usage and deletions

### `geneUsage()`

```r
geneUsage(clip_db, chain = c("IGH", "IGK", "IGL", "TRB"),
          genes_order = NULL, rmPseudo = TRUE)
```

Per-gene usage frequencies, used as input to the deletion tests.

### `deletionsByBinom()`

Double-chromosome deletion detection across a population via a binomial test on gene usage.

### `deletionsByVpooled()`

Single-chromosome deletion detection for V genes using a pooled approach.

### `nonReliableVGenes()`

Identify V genes that are unreliable (frequently appearing in low-confidence multi-assignments), e.g.
in partial-coverage data. The result is passed to `deletionsByBinom()` and `createFullHaplotype()`.

## Visualization

| Function | Output |
| --- | --- |
| `plotHaplotype()` | Per-subject haplotype map (static, or interactive with `html_output = TRUE`) |
| `hapHeatmap()` | Multi-sample allele-assignment heatmap; returns `list(p, width, height)` |
| `hapDendo()` | Dendrogram + heatmap clustering samples by haplotype similarity (needs `ggdendro`) |
| `deletionHeatmap()` | Single-chromosome deletion heatmap across subjects (static or interactive) |
| `plotDeletionsByBinom()` | Double-chromosome deletions across the population |
| `plotDeletionsByVpooled()` | V/D/J single-chromosome deletions |

Common arguments: `chain`, `genes_order` (override gene ordering), `removeIGH`/`removeIGH`-style label
trimming, `lk_cutoff`/`kThreshDel` (confidence thresholds), and `html_output` for the interactive
variants (requires `plotly` + `htmlwidgets`).

## Bundled data

`samples_db` (example IGH repertoire), `samplesHaplotype` (example result), the germline references
`HVGERM`/`HDGERM`/`HJGERM`/`KVGERM`/`KJGERM`/`LVGERM`/`LJGERM`, and `GENE.loc` (chromosomal gene
ordering). See [Input format](input.md).
