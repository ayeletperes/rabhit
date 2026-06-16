# Input format

RAbHIT works on pre-processed antibody (or TRB) repertoire sequencing data with heterozygosity in at
least one gene, plus a database of germline gene sequences.

## Repertoire data

A `data.frame` in [AIRR format](https://docs.airr-community.org/), where each row is a unique
observation. The required columns are:

| Column | Description |
| --- | --- |
| `subject` | Subject identifier |
| `v_call` | (Comma-separated) name(s) of the nearest V allele(s), IMGT format |
| `d_call` | (Comma-separated) name(s) of the nearest D allele(s), IMGT format |
| `j_call` | (Comma-separated) name(s) of the nearest J allele(s), IMGT format |

!!! note
    Allele names are expected in IMGT nomenclature (`GENE*ALLELE`, e.g. `IGHV1-2*01`). Support for
    PIgLET allele-similarity-cluster (ASC) names and OGRDB reference sets is provided via integration
    with the [`piglet`](https://cran.r-project.org/package=piglet) package.

## Germline database

A named character vector of IMGT-gapped germline nucleotide sequences, where names are the IMGT allele
calls (e.g. `IGHV1-2*01`). RAbHIT bundles reference databases for the supported chains:

| Object | Segment |
| --- | --- |
| `HVGERM`, `HDGERM`, `HJGERM` | IGH V / D / J |
| `KVGERM`, `KJGERM` | IGK V / J |
| `LVGERM`, `LJGERM` | IGL V / J |

Load them with `data(HVGERM, HDGERM, HJGERM)`. If you do not pass `toHap_GERM`, RAbHIT selects the
appropriate bundled germline for the requested `chain`.

## Example data

The package ships an example repertoire (`samples_db`) and a pre-computed haplotype result
(`samplesHaplotype`) used throughout the [Tutorial](tutorial.md):

```r
data(samples_db, samplesHaplotype)
```
