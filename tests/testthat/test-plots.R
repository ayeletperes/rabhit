# Smoke tests for the plotting functions on the bundled haplotype result.
# The interactive (plotly) paths depend on Suggested packages and are skipped
# when those are not installed, which also exercises the requireNamespace guards.

test_that("plotHaplotype produces a static plot", {
  data(samplesHaplotype, envir = environment())
  hap <- samplesHaplotype[samplesHaplotype$subject == "I5", ]
  expect_error(plotHaplotype(hap, html_output = FALSE), NA)
})

test_that("plotHaplotype(html_output = TRUE) needs plotly, else errors gracefully", {
  data(samplesHaplotype, envir = environment())
  hap <- samplesHaplotype[samplesHaplotype$subject == "I5", ]
  if (requireNamespace("plotly", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    p <- plotHaplotype(hap, html_output = TRUE)
    expect_true(inherits(p, "plotly") || inherits(p, "htmlwidget"))
  } else {
    expect_error(plotHaplotype(hap, html_output = TRUE), "plotly")
  }
})

test_that("hapHeatmap returns a plot with sizing metadata", {
  data(samplesHaplotype, envir = environment())
  res <- hapHeatmap(samplesHaplotype)
  expect_type(res, "list")
  expect_true(all(c("width", "height") %in% names(res)))
})

test_that("deletionHeatmap produces a static plot", {
  data(samplesHaplotype, envir = environment())
  expect_error(deletionHeatmap(samplesHaplotype, html_output = FALSE), NA)
})

test_that("hapDendo needs ggdendro, else errors gracefully", {
  data(samplesHaplotype, envir = environment())
  if (requireNamespace("ggdendro", quietly = TRUE)) {
    expect_error(hapDendo(samplesHaplotype), NA)
  } else {
    expect_error(hapDendo(samplesHaplotype), "ggdendro")
  }
})
