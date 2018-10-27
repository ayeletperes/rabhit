load(file.path('..','data-tests','TestDb.rda'))

test_that("Test creatFullHaplotype",{

  # compare number of deletion inferred
  hap_df <- createFullHaplotype(sample_db,toHap_GERM = c(VGERM,DGERM), supress_print = T)
  del <- nrow(hap_df[grepl('Del',hap_df$IGHJ6_02) | grepl('Del',hap_df$IGHJ6_03),])
  expect_equal(del, 13)

  # check that the germlime matches the genes to haplotype
  expect_error(createFullHaplotype(sample_db,toHap_GERM = c(IGKV), supress_print = T),"Genes in haplotype column to be infered do not match the genes germline given")
})

