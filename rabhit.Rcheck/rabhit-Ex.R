pkgname <- "rabhit"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "rabhit-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('rabhit')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("createFullHaplotype")
### * createFullHaplotype

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: createFullHaplotype
### Title: Anchor gene haplotype inference
### Aliases: createFullHaplotype

### ** Examples

# Load example data and germlines
data(samples_db, HVGERM, HDGERM)

# Selecting a single individual
clip_db = samples_db[samples_db$subject=='I5', ]

# Infering haplotype
haplo_db = createFullHaplotype(clip_db,toHap_col=c('v_call','d_call'),
hapBy_col='j_call',hapBy='IGHJ6',toHap_GERM=c(HVGERM,HDGERM))





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("createFullHaplotype", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("deletionHeatmap")
### * deletionHeatmap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: deletionHeatmap
### Title: Graphical output of single chromosome deletions
### Aliases: deletionHeatmap

### ** Examples

# Plotting single choromosme deletion from haplotype inference
deletionHeatmap(samplesHaplotype)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("deletionHeatmap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("deletionsByBinom")
### * deletionsByBinom

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: deletionsByBinom
### Title: Double chromosome deletion by relative gene usage
### Aliases: deletionsByBinom

### ** Examples

# Load example data and germlines
data(samples_db)

# Selecting a single individual
clip_db = samples_db[samples_db$subject=='I5', ]
# Infering haplotype
del_binom_df = deletionsByBinom(clip_db)
head(del_binom_df)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("deletionsByBinom", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("deletionsByVpooled")
### * deletionsByVpooled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: deletionsByVpooled
### Title: Single chromosomal D or J gene deletions inferred by the V
###   pooled method
### Aliases: deletionsByVpooled

### ** Examples

## No test: 
data(samples_db)

# Infering V pooled deletions
del_db <- deletionsByVpooled(samples_db)
head(del_db)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("deletionsByVpooled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hapDendo")
### * hapDendo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hapDendo
### Title: Hierarchical clustering of haplotypes graphical output
### Aliases: hapDendo

### ** Examples

# Plotting haplotype hierarchical clustering based on the Jaccard distance
## No test: 
hapDendo(samplesHaplotype)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hapDendo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hapHeatmap")
### * hapHeatmap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hapHeatmap
### Title: Graphical output of alleles division by chromosome
### Aliases: hapHeatmap

### ** Examples

# Plotting haplotpe heatmap
p <- hapHeatmap(samplesHaplotype)
p$p



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hapHeatmap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("nonReliableVGenes")
### * nonReliableVGenes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: nonReliableVGenes
### Title: Detect non reliable gene assignment
### Aliases: nonReliableVGenes

### ** Examples

# Example IGHV call data frame
clip_db <- data.frame(subject=rep('S1',6),
v_call=c('IGHV1-69*01','IGHV1-69*01','IGHV1-69*01,IGHV1-69*02',
'IGHV4-59*01,IGHV4-61*01','IGHV4-59*01,IGHV4-31*02','IGHV4-59*01'))
# Detect non reliable genes
nonReliableVGenes(clip_db)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("nonReliableVGenes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotDeletionsByBinom")
### * plotDeletionsByBinom

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotDeletionsByBinom
### Title: Graphical output of double chromosome deletions
### Aliases: plotDeletionsByBinom

### ** Examples


# Load example data and germlines
data(samples_db)

# Infering haplotype
deletions_db = deletionsByBinom(samples_db);
plotDeletionsByBinom(deletions_db)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotDeletionsByBinom", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotDeletionsByVpooled")
### * plotDeletionsByVpooled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotDeletionsByVpooled
### Title: Graphical output for single chromosome D or J gene deletions
###   according to V pooled method
### Aliases: plotDeletionsByVpooled

### ** Examples

## No test: 
# Load example data and germlines
data(samples_db)
del_db <- deletionsByVpooled(samples_db)
plotDeletionsByVpooled(del_db)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotDeletionsByVpooled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotHaplotype")
### * plotHaplotype

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotHaplotype
### Title: Graphical output of an inferred haplotype
### Aliases: plotHaplotype

### ** Examples


# Selecting a single individual from the haplotype samples data
haplo_db = samplesHaplotype[samplesHaplotype$subject=='I5', ]

# plot haplotype
plotHaplotype(haplo_db)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotHaplotype", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
