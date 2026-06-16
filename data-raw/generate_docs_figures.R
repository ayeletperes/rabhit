# Generate the figures embedded in the MkDocs documentation site.
# The docs workflow (docs.yaml) is Python-only and cannot run R, so the example
# plots are rendered to PNG here and committed under site_docs/assets/img/.
#
# Run from the package root with:  Rscript data-raw/generate_docs_figures.R

library(rabhit)

out_dir <- "site_docs/assets/img"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data(samples_db, HVGERM, HDGERM)

png_file <- function(name) file.path(out_dir, name)

## 1. Single-subject haplotype map (IGHJ6 anchor) -----------------------------
clip_db <- samples_db[samples_db$subject == "I5", ]
haplo_J6 <- createFullHaplotype(clip_db,
                                toHap_col = c("v_call", "d_call"),
                                hapBy_col = "j_call", hapBy = "IGHJ6",
                                toHap_GERM = c(HVGERM, HDGERM))

png(png_file("haplotype_single.png"), width = 1500, height = 1100, res = 120)
plotHaplotype(haplo_J6)
dev.off()

## Multi-subject haplotype table (reused below) -------------------------------
clip_dbs <- samples_db[samples_db$subject != "I5_FR2", ]
haplo_multi <- createFullHaplotype(clip_dbs,
                                   toHap_col = c("v_call", "d_call"),
                                   hapBy_col = "j_call", hapBy = "IGHJ6",
                                   toHap_GERM = c(HVGERM, HDGERM))

## 2. Haplotype heatmap across subjects ---------------------------------------
p.list <- hapHeatmap(haplo_multi)
png(png_file("hap_heatmap.png"),
    width = max(1400, p.list$width * 80),
    height = max(1000, p.list$height * 80), res = 120)
print(p.list$p)
dev.off()

## 3. Haplotype dendrogram ----------------------------------------------------
png(png_file("hap_dendro.png"), width = 1600, height = 900, res = 120)
print(hapDendo(haplo_multi))
dev.off()

## 4. Deletion heatmap (population) -------------------------------------------
nonReliable_Vgenes <- nonReliableVGenes(samples_db)
del_binom_db <- deletionsByBinom(samples_db, nonReliable_Vgenes = nonReliable_Vgenes)
haplo_del <- createFullHaplotype(samples_db,
                                 toHap_col = c("v_call", "d_call"),
                                 hapBy_col = "j_call", hapBy = "IGHJ6",
                                 toHap_GERM = c(HVGERM, HDGERM),
                                 deleted_genes = del_binom_db,
                                 nonReliable_Vgenes = nonReliable_Vgenes)
png(png_file("deletion_heatmap.png"), width = 1600, height = 1000, res = 120)
print(deletionHeatmap(haplo_del))
dev.off()

## 5. Double-chromosome deletions (binomial) ----------------------------------
del_plot_db <- del_binom_db[grep("IGHJ", del_binom_db$gene, invert = TRUE), ]
png(png_file("deletions_binom.png"), width = 1600, height = 950, res = 120)
print(plotDeletionsByBinom(del_plot_db))
dev.off()

## 6. Single-chromosome V deletions (pooled) ----------------------------------
del_db <- deletionsByVpooled(samples_db, nonReliable_Vgenes = nonReliable_Vgenes)
png(png_file("deletions_vpooled.png"), width = 1400, height = 700, res = 120)
print(plotDeletionsByVpooled(del_db))
dev.off()

message("Figures written to ", out_dir)
