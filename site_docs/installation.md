# Installation

## From CRAN

Once RAbHIT is available on CRAN:

```r
install.packages("rabhit")
```

## Development version (GitHub)

```r
# install.packages("devtools")
devtools::install_github("ayeletperes/rabhit")
```

## From source

Install the build dependencies, then clone and build:

```r
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "plotly"))
```

```r
library(devtools)
install_deps()
document()
build()
install()
```

## Optional features

Some functionality lives behind optional (Suggested) packages so the core install stays light. Install
them only if you need the corresponding feature:

| Feature | Package(s) |
| --- | --- |
| Interactive HTML5 plots (`html_output = TRUE`) | `plotly`, `htmlwidgets` |
| Haplotype dendrogram (`hapDendo`) | `ggdendro` |
| ASC / OGRDB allele-cluster naming | `piglet` |

If a function that needs one of these is called without it installed, RAbHIT stops with an informative
message telling you what to install.
