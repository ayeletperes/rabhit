# Allele clusters (ASC / OGRDB)

Some IGHV alleles are nearly identical and are shared between duplicated genes
(e.g. *IGHV1-69* / *IGHV1-69D*, *IGHV3-23* / *IGHV3-23D*). These cannot always be
assigned unambiguously from short reads, which weakens genotype and haplotype
inference. **Allele Similarity Clusters (ASCs)**, from the
[PIgLET](https://cran.r-project.org/package=piglet) package, group such alleles
into clusters with a population-informed threshold and give them a consistent
naming scheme (`IGHVF<family>-G<cluster>`) that is compatible with
[OGRDB](https://ogrdb.airr-community.org) reference sets.

RAbHIT integrates this through a single helper, [`convertToASC()`](reference.md),
which prepares ASC-named data and germline for the standard haplotype workflow.

!!! note "Optional dependency"
    This feature requires the `piglet` package (and a sequence-distance backend such
    as `DECIPHER`). Install it with `install.packages("piglet")`. RAbHIT only loads it
    when you call `convertToASC()`.

## Basic usage

```r
library(rabhit)
data(samples_db, HVGERM)

clip_db <- samples_db[samples_db$subject == "I5", ]

# Convert V calls + germline to ASC names (clusters inferred from the germline)
asc <- convertToASC(clip_db, HVGERM, chain = "IGH")

# Haplotype the ASC-named V genes, anchored on IGHJ6
hap <- createFullHaplotype(asc$clip_db, toHap_col = "v_call",
                           hapBy_col = "j_call", hapBy = "IGHJ6",
                           toHap_GERM = asc$germline, chain = "IGH")

# Plot with the ASC gene order returned by convertToASC()
plotHaplotype(hap, genes_order = asc$genes_order)
```

`convertToASC()` returns a list with the ASC-named `clip_db`, the ASC `germline`
(pass it as `toHap_GERM`), a `genes_order` for the plotting functions, and the
`allele_cluster_table` that was used.

## Reproducible / OGRDB clusters

By default the cluster table is **inferred from the germline**, which can vary with
the reference set. For reproducible results — and to align with the community
reference — supply a fixed ASC table, e.g. the published thresholds from OGRDB/Zenodo:

```r
# Download the published ASC archive and extract its threshold table
archive <- piglet::recentAlleleClusters(get_file = TRUE)
asc_table <- piglet::extractASCTable(archive)

asc <- convertToASC(clip_db, HVGERM, chain = "IGH",
                    allele_cluster_table = asc_table)
```

You can also pass a table produced once with `piglet::inferAlleleClusters()` and
reuse it across subjects.

## Notes

- Clustering is applied to **V genes only**; D and J calls are left unchanged, so
  D/J haplotyping continues to use standard IMGT names.
- Duplicated genes are largely resolved for free: alleles shared across duplicated
  genes fall into the same cluster, removing the cross-assignment ambiguity that
  RAbHIT would otherwise flag as non-reliable.
