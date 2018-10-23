## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
library(rabhit)
# Load example sequence data and example germline database
data(sample_db, HVGERM, HDGERM)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Infered haplotype summary table
haplo_db <- createFullHaplotype(sample_db,toHap_col=c("V_CALL","D_CALL"),
hapBy_col="J_CALL",hapBy="IGHJ6",toHap_GERM=c(HVGERM,HDGERM))

head(haplo_db)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Plot the haplotype
plotHaplotype(haplo_db)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Plot interactive haplotype plot
p <- plotHaplotype(haplo_db,html_output = T)
#save plot to html output
htmlwidgets::saveWidget(p, "haplotype.html",selfcontained = T)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Infered deletion summary table
del_binom_db <- deletionsByBinom(samples_db)
head(del_binom_db)

## ---- eval=TRUE, warning=FALSE,fig.height=9,fig.width=12-----------------
# Infered deletion summary table
plotDeletionsByBinom(del_binom_db)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Infered deletion summary table
del_binom_db <- deletionsByBinom(samples_db)
haplo_db <- createFullHaplotype(sample_db,toHap_col=c("V_CALL","D_CALL"),
hapBy_col="J_CALL",hapBy="IGHJ6",toHap_GERM=c(HVGERM,HDGERM),deleted_genes = del_binom_db,supress_print = T)
plotHaplotype(haplo_db)

## ---- eval=TRUE, warning=FALSE,fig.height=16,fig.width=15----------------
# Load example sequence data
data(samples_db)
# Infered haplotype summary table for multiple subjects
haplo_db <- createFullHaplotype(samples_db,toHap_col=c("V_CALL","D_CALL"),
hapBy_col="J_CALL",hapBy="IGHJ6",toHap_GERM=c(HVGERM,HDGERM))
# plot deletion heatmap
deletionHeatmap(haplo_db)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Infered deletion summary table
del_db <- deletionsByVpooled(samples_db)
head(del_db)

## ---- eval=TRUE, warning=FALSE,fig.height=4,fig.width=8------------------
# Plot the deletion heatmap
plotDeletionsByVpooled(del_db)

