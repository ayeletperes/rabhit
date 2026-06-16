# RAbHIT — R Antibody Haplotype Inference Tool

Analysis of antibody repertoires by high-throughput sequencing is of major importance in
understanding adaptive immune responses. Our knowledge of variations in the genomic loci encoding
antibody genes is incomplete, mostly due to the technical difficulty of aligning short reads to
these highly repetitive loci. This partial knowledge results in conflicting V-D-J gene assignments
between different algorithms, and biased genotype and haplotype inference. Previous studies have
shown that haplotypes can be inferred by taking advantage of IGHJ6 heterozygosity, observed in
approximately one third of the population.

**RAbHIT is a haplotype inference tool based on a robust method for determining V-D-J haplotypes by
adapting a Bayesian framework.** The method extends haplotype inference to IGHD, IGHV, IGKJ, IGKV,
and IGLV based analysis, thereby enabling inference of complex genetic events like deletions and
copy-number variations across the population. It also supports the TRB chain. For each haplotyped
gene the tool reports a **Bayes factor** (`lK = log10` of the likelihood ratio) that quantifies the
certainty of the inference.

## Core abilities

- **Haplotype inference** from a J, D, or V gene used as anchor
- **Single-chromosome deletion** detection
- **Double-chromosome deletion** detection
- Publication-quality static and interactive (HTML5) **visualizations**

## Supported chains

`IGH`, `IGK`, `IGL`, and `TRB`.

## At a glance

```r
library(rabhit)
data(samples_db, HVGERM, HDGERM)

# one individual, haplotype by IGHJ6
clip_db <- samples_db[samples_db$subject == "I5", ]
haplo <- createFullHaplotype(clip_db,
                             toHap_col = c("v_call", "d_call"),
                             hapBy_col = "j_call", hapBy = "IGHJ6",
                             toHap_GERM = c(HVGERM, HDGERM))
plotHaplotype(haplo)
```

See the [Tutorial](tutorial.md) for the full workflow and the
[Function reference](reference.md) for the complete API.

## Citation

If you use RAbHIT, please cite the papers listed on the [Citation](citation.md) page.

RAbHIT is free for use under the [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
license.
