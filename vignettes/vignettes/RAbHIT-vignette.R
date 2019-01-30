## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
library(rabhit)
# Load example sequence data and example germline database
data(sample_db, HVGERM, HDGERM)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Inferred haplotype summary table
haplo_db <- createFullHaplotype(sample_db, toHap_col=c("V_CALL","D_CALL"),
                                hapBy_col="J_CALL", hapBy="IGHJ6", 
                                toHap_GERM=c(HVGERM, HDGERM))

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  head(haplo_db)

## ---- eval=TRUE, warning=FALSE,echo=FALSE--------------------------------
head(haplo_db,3)

## ---- eval=TRUE, warning=FALSE,fig.width=15,fig.height=10----------------
# Plot the haplotype map
plotHaplotype(haplo_db)

## ---- eval=TRUE, warning=FALSE, fig.width=12, fig.height=8---------------
# Load example sequence data
data(samples_db)
# Inferred haplotype summary table
haplo_db <- createFullHaplotype(samples_db, toHap_col=c("V_CALL","D_CALL"),
                        hapBy_col="J_CALL", hapBy="IGHJ6", toHap_GERM=c(HVGERM, HDGERM),
                        supress_print = T)
# Plot the haplotype inferred heatmap
hapHeatmap(haplo_db)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Inferred deletion summary table
del_binom_db <- deletionsByBinom(samples_db)
head(del_binom_db)

## ---- eval=TRUE, warning=FALSE,fig.height=9,fig.width=15-----------------
# Don't plot IGHJ
del_binom_db <- del_binom_db[grep('IGHJ',del_binom_db$GENE, invert = T),]
# Inferred deletion summary table
plotDeletionsByBinom(del_binom_db) 

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Inferred deletion summary table
del_binom_db <- deletionsByBinom(sample_db)
# Input the summary table to the deleted_genes flag
haplo_db <- createFullHaplotype(sample_db, toHap_col=c("V_CALL","D_CALL"),
                        hapBy_col="J_CALL", hapBy="IGHJ6", toHap_GERM=c(HVGERM, HDGERM),
                        deleted_genes = del_binom_db, supress_print = T)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Generate interactive haplotype plot
p <- plotHaplotype(haplo_db, html_output = TRUE)

# Save plot to html output
htmlwidgets::saveWidget(p, "haplotype.html", selfcontained = T)

# Plot the haplotype
p

## ---- eval=TRUE, warning=FALSE,fig.height=12,fig.width=15----------------
# Inferred haplotype summary table for multiple subjects
haplo_db <- createFullHaplotype(samples_db, toHap_col=c("V_CALL","D_CALL"), 
                                hapBy_col="J_CALL", hapBy="IGHJ6", 
                                toHap_GERM=c(HVGERM, HDGERM), supress_print = T)
# plot deletion heatmap
deletionHeatmap(haplo_db)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Inferred deletion summary table
del_db <- deletionsByVpooled(samples_db)
head(del_db)

## ---- eval=TRUE, warning=FALSE,fig.height=4,fig.width=8------------------
# Plot the deletion heatmap
plotDeletionsByVpooled(del_db)

