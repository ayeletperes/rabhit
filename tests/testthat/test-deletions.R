# Tests for gene-usage and deletion-detection helpers.

test_that("geneUsage returns per-gene usage frequencies", {
  data(samples_db, envir = environment())
  clip_db <- samples_db[samples_db$subject == "I5", ]
  gu <- geneUsage(clip_db, chain = "IGH")
  expect_s3_class(gu, "data.frame")
  expect_true("gene" %in% names(gu))
  expect_gt(nrow(gu), 0)
})

test_that("nonReliableVGenes returns gene names for a partial-coverage sample", {
  data(samples_db, envir = environment())
  clip_db <- samples_db[samples_db$subject == "I5_FR2", ]
  nr <- nonReliableVGenes(clip_db)
  # Returns a (possibly empty) list/character of non-reliable V genes.
  expect_true(is.list(nr) || is.character(nr) || is.null(nr))
})

test_that("deletionsByBinom detects double-chromosome deletions", {
  data(samples_db, envir = environment())
  clip_db <- samples_db[samples_db$subject == "I5", ]
  del <- deletionsByBinom(clip_db, chain = "IGH")
  expect_s3_class(del, "data.frame")
  expect_true("gene" %in% names(del))
  expect_gt(nrow(del), 0)
})

test_that("deletionsByVpooled runs on a multi-subject set", {
  skip_on_cran()
  data(samples_db, envir = environment())
  clip_db <- samples_db[samples_db$subject %in% c("I1", "I2"), ]
  del <- deletionsByVpooled(clip_db, chain = "IGH")
  expect_s3_class(del, "data.frame")
  expect_true("gene" %in% names(del))
})
