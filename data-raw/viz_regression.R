# Visualization regression harness.
# Renders every plot function from the bundled data to PNGs (and exercises the
# plotly/HTML branches) so output can be compared before/after the refactor.
#
# Usage:  Rscript data-raw/viz_regression.R <out_dir>
# Loads the *source* tree (devtools::load_all) so it reflects working changes.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "/tmp/rabhit_viz/current"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(devtools::load_all(".", quiet = TRUE))

data(samples_db, HVGERM, HDGERM)
png_file <- function(name) file.path(out_dir, name)

## 1. Single-subject haplotype map (IGHJ6 anchor) -----------------------------
clip_db <- samples_db[samples_db$subject == "I5", ]
haplo_J6 <- createFullHaplotype(clip_db,
                                toHap_col = c("v_call", "d_call"),
                                hapBy_col = "j_call", hapBy = "IGHJ6",
                                toHap_GERM = c(HVGERM, HDGERM))
png(png_file("01_haplotype_single.png"), width = 1500, height = 1100, res = 120)
plotHaplotype(haplo_J6)
dev.off()

## plotly/HTML branch of plotHaplotype -----------------------------------------
hp <- plotHaplotype(haplo_J6, html_output = TRUE)
saveRDS(class(hp), png_file("01b_haplotype_html_class.rds"))

## Multi-subject haplotype table ----------------------------------------------
clip_dbs <- samples_db[samples_db$subject != "I5_FR2", ]
haplo_multi <- createFullHaplotype(clip_dbs,
                                   toHap_col = c("v_call", "d_call"),
                                   hapBy_col = "j_call", hapBy = "IGHJ6",
                                   toHap_GERM = c(HVGERM, HDGERM))

## 2. Haplotype heatmap -------------------------------------------------------
p.list <- hapHeatmap(haplo_multi)
png(png_file("02_hap_heatmap.png"),
    width = max(1400, p.list$width * 80),
    height = max(1000, p.list$height * 80), res = 120)
print(p.list$p)
dev.off()

## 3. Haplotype dendrogram ----------------------------------------------------
png(png_file("03_hap_dendro.png"), width = 1600, height = 900, res = 120)
print(hapDendo(haplo_multi))
dev.off()

## 4. Deletion heatmap --------------------------------------------------------
nonReliable_Vgenes <- nonReliableVGenes(samples_db)
del_binom_db <- deletionsByBinom(samples_db, nonReliable_Vgenes = nonReliable_Vgenes)
haplo_del <- createFullHaplotype(samples_db,
                                 toHap_col = c("v_call", "d_call"),
                                 hapBy_col = "j_call", hapBy = "IGHJ6",
                                 toHap_GERM = c(HVGERM, HDGERM),
                                 deleted_genes = del_binom_db,
                                 nonReliable_Vgenes = nonReliable_Vgenes)
png(png_file("04_deletion_heatmap.png"), width = 1600, height = 1000, res = 120)
print(deletionHeatmap(haplo_del))
dev.off()

## 5. Double-chromosome deletions (binomial) ----------------------------------
del_plot_db <- del_binom_db[grep("IGHJ", del_binom_db$gene, invert = TRUE), ]
png(png_file("05_deletions_binom.png"), width = 1600, height = 950, res = 120)
print(plotDeletionsByBinom(del_plot_db))
dev.off()

## 6. Single-chromosome V deletions (pooled) ----------------------------------
del_db <- deletionsByVpooled(samples_db, nonReliable_Vgenes = nonReliable_Vgenes)
png(png_file("06_deletions_vpooled.png"), width = 1400, height = 700, res = 120)
print(plotDeletionsByVpooled(del_db))
dev.off()

message("Wrote regression renders to ", out_dir)
