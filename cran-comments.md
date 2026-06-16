## Resubmission

This is a resubmission of 'rabhit', which was archived on CRAN on 2024-03-26.
In this version we have:

* Reduced the dependency footprint to lower the maintenance burden and the
  risk that an upstream change again breaks the package:
  - Removed the hard dependencies on 'gtools', 'RColorBrewer', 'dendextend',
    'rlang' and 'methods'. The two small pieces that were used (a Dirichlet
    density and a ColorBrewer palette) are now inlined; 'dendextend', 'rlang'
    and 'methods' were unused imports.
  - Moved 'plotly', 'htmlwidgets' and 'ggdendro' to Suggests; the functions
    that need them now check for them with `requireNamespace()` and fail
    gracefully with an informative message when they are not installed.
  - Moved 'ggplot2' from Depends to Imports and declared the base/recommended
    packages it imports from ('stats', 'utils', 'graphics', 'grDevices') in
    Imports.
* Removed stray build/check artifacts that had been committed to the source
  tree (e.g. a top-level `Rplots.pdf`).
* Added a 'testthat' test suite.

## Test environments

* local: Ubuntu, R 4.5.1
* (to be completed) win-builder (devel and release)
* (to be completed) macOS / R-hub

## R CMD check results

0 errors | 0 warnings | 0-1 notes

The only NOTE seen locally is the environment-specific
"unable to verify current time", which is not related to the package.

## Downstream dependencies

There are currently no reverse dependencies for this package.
