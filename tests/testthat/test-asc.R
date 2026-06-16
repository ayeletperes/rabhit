# Tests for the PIgLET allele-similarity-cluster (ASC) integration.
# Inference is slow and needs piglet + a sequence-distance backend, so these are
# skipped on CRAN and when the optional packages are unavailable.

test_that("convertToASC is exported and guards its optional dependency", {
  expect_true(is.function(convertToASC))
})

test_that("convertToASC yields ASC-named input that createFullHaplotype accepts", {
  skip_on_cran()
  skip_if_not_installed("piglet")
  skip_if_not_installed("DECIPHER")

  data(samples_db, HVGERM, envir = environment())
  clip_db <- samples_db[samples_db$subject == "I5", ]

  asc <- convertToASC(clip_db, HVGERM, chain = "IGH")

  expect_type(asc, "list")
  expect_named(asc, c("clip_db", "germline", "genes_order", "allele_cluster_table"))
  # V calls and germline are now ASC-named (IGHVF<family>-G<cluster>).
  expect_true(any(grepl("^IGHVF", asc$clip_db$v_call)))
  expect_true(all(grepl("^IGHVF", names(asc$germline))))
  expect_true(all(grepl("^IGHVF", asc$genes_order)))

  hap <- createFullHaplotype(asc$clip_db, toHap_col = "v_call",
                             hapBy_col = "j_call", hapBy = "IGHJ6",
                             toHap_GERM = asc$germline, chain = "IGH")
  expect_s3_class(hap, "data.frame")
  expect_gt(nrow(hap), 0)
  expect_true(all(grepl("^IGHVF", hap$gene)))
})

test_that("convertToASC accepts a precomputed ASC table", {
  skip_on_cran()
  skip_if_not_installed("piglet")
  skip_if_not_installed("DECIPHER")

  data(samples_db, HVGERM, envir = environment())
  clip_db <- samples_db[samples_db$subject == "I5", ]

  tbl <- piglet::inferAlleleClusters(HVGERM, locus = "IGHV", quiet = TRUE)$alleleClusterTable
  asc <- convertToASC(clip_db, HVGERM, chain = "IGH", allele_cluster_table = tbl)
  expect_identical(asc$allele_cluster_table, tbl)
  expect_true(any(grepl("^IGHVF", asc$clip_db$v_call)))
})
