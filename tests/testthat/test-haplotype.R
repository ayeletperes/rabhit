# Tests for the core haplotype inference (createFullHaplotype / createHaplotypeTable).
# These also act as a regression guard: the result for subject I5 must match the
# bundled `samplesHaplotype` reference, which protects the inlined Dirichlet density
# (replacing gtools::ddirichlet) against numerical drift.

test_that("createFullHaplotype reproduces the bundled reference (J6 anchor)", {
  data(samples_db, HVGERM, HDGERM, samplesHaplotype, envir = environment())
  clip_db <- samples_db[samples_db$subject == "I5", ]
  hap <- createFullHaplotype(clip_db, toHap_col = c("v_call", "d_call"),
                             hapBy_col = "j_call", hapBy = "IGHJ6",
                             toHap_GERM = c(HVGERM, HDGERM))

  expect_s3_class(hap, "data.frame")
  expect_true(all(c("subject", "gene", "IGHJ6_02", "IGHJ6_03", "alleles") %in% names(hap)))
  expect_gt(nrow(hap), 0)

  ref <- samplesHaplotype[samplesHaplotype$subject == "I5", ]
  m <- merge(ref, hap, by = c("subject", "gene"), suffixes = c(".ref", ".new"))
  expect_gt(nrow(m), 0)
  expect_identical(m$IGHJ6_02.new, m$IGHJ6_02.ref)
  expect_identical(m$IGHJ6_03.new, m$IGHJ6_03.ref)
})

test_that("createFullHaplotype works with a D gene as anchor", {
  data(samples_db, HVGERM, envir = environment())
  clip_db <- samples_db[samples_db$subject == "I5", ]
  hap <- createFullHaplotype(clip_db, toHap_col = "v_call",
                             hapBy_col = "d_call", hapBy = "IGHD2-21",
                             toHap_GERM = HVGERM)

  expect_s3_class(hap, "data.frame")
  expect_true(any(grepl("^IGHV", hap$gene)))
})

test_that("createFullHaplotype skips an anchor that is not a 2-allele heterozygote", {
  data(samples_db, HVGERM, HDGERM, envir = environment())
  clip_db <- samples_db[samples_db$subject == "I5", ]
  # A gene that is not present in the data cannot yield two anchor alleles, so the
  # sample is reported and skipped rather than haplotyped.
  expect_message(
    res <- createFullHaplotype(clip_db, toHap_col = c("v_call", "d_call"),
                               hapBy_col = "j_call", hapBy = "IGHJ999",
                               toHap_GERM = c(HVGERM, HDGERM)),
    "two alleles"
  )
  expect_true(is.null(res) || nrow(res) == 0)
})
