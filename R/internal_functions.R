# Internal functions -----------------------------------------------------

#' @include rabhit.R
NULL

# Calculate models likelihood
#
# \code{createHaplotypeTable} calculate likelihoods
#
# @param    X      a vector of counts
# @param    alpha_dirichlet      alpha parameter for dirichlet distribution
# @param    epsilon    epsilon
# @param    priors      a vector of priors
#
# @return  log10 of the likelihoods
#
#
get_probabilites_with_priors <- function(X,alpha_dirichlet=c(0.5,0.5)*2,epsilon=0.01,params=c(0.5,0.5)){
  ## Hypotheses
  X<-sort(X,decreasing=TRUE)
  Number_Of_Divisions<-0

  H1<-c(1,0)
  H2<-c(params[1],params[2])

  E1<-ddirichlet((H1+epsilon)/sum(H1+epsilon),alpha_dirichlet+X)
  E2<-ddirichlet((H2+epsilon)/sum(H2+epsilon),alpha_dirichlet+X)

  while(sort(c(E1,E2),decreasing=TRUE)[2] == 0 ){
    Number_Of_Divisions<-Number_Of_Divisions+1
    X <- X/10
    E1<-ddirichlet((H1+epsilon)/sum(H1+epsilon),alpha_dirichlet+X)
    E2<-ddirichlet((H2+epsilon)/sum(H2+epsilon),alpha_dirichlet+X)
  }
  return(c(log10(c(E1,E2)),Number_Of_Divisions))
}

##############################################################################################################

# Create haplotype table
#
# \code{createHaplotypeTable} Haplotype of a specific gene
#
# @details
#
# @param  df  table of counts
# @param  HapByPriors vector of frequencies of each of the anchor gene alleles
# @param  toHapByCol logical, haplotype each chromosome separetly to imrove the aligner assignmnet
# @param  toHapPriors vector of frequencies of the haplotyped gene alleles
#
# @return  data frame with chromosomal associasions of alleles of a specific gene
#
# @examples
# # Load example data and germlines
# data(sample_db)
# data(germline_igh)
#
#
# @export
createHaplotypeTable <- function(df,HapByPriors=c(0.5,0.5),toHapByCol=TRUE,toHapPriors=c(0.5,0.5)){
  hapBy <- colnames(df)
  tohap <- rownames(df)
  tohap.gene <- strsplit(tohap[1],'*',fixed = T)[[1]][1]

  GENES.df <- data.frame(GENE=tohap.gene,'Unk','Unk',stringsAsFactors = F)
  names(GENES.df)[2:3] <- gsub('*','.',hapBy,fixed = T)
  GENES.df.num <- melt(df)

  df.old <- df
  if(toHapByCol){
    if(nrow(df) > 1){
      for(j in 1:ncol(df)){

        counts <- c(sort(df[,j],decreasing = T),0,0,0)
        if(sum(counts) != counts[1]){
          names(counts)[1:nrow(df)] <- names(sort(df[,j],decreasing = T))

          toHapPriors_srtd <- if(!is.null(names(toHapPriors))) toHapPriors[names(counts[1:2])] else toHapPriors
          resCol <- get_probabilites_with_priors(counts[1:2],params = toHapPriors_srtd)
          resMaxInd <- which.max(resCol[-(length(resCol))])
          if(resMaxInd < nrow(df)){
            df[,j][order( df[,j],decreasing = T)[(resMaxInd+1):nrow(df)]] <- 0
          }
        }

      }
    }

  }

  counts.list <- list()
  res.list <- list()

  for(i in 1:nrow(df)){
    allele <- rownames(df)[i]
    gene <- strsplit(allele,'*',fixed=T)[[1]][1]


    counts <- c(sort(df[i,],decreasing = T),0,0,0)

    if(ncol(df)==1) names(counts)[1:ncol(df)] <- colnames(df )
    if(ncol(df)!=1) names(counts)[1:ncol(df)] <- names(sort(df[i,],decreasing = T))

    HapByPriors_srtd <- if(!is.null(names(HapByPriors))) HapByPriors[names(counts[1:2])] else HapByPriors
    res <- get_probabilites_with_priors(counts[1:2],params = HapByPriors_srtd)


    # Assign allele for a chromosome
    ### If hetero in chomosome : check if anythong was assigned (equals 'unk'), if was (different than 'unk'), paste second allele in the same chromosome
    if(res[1] > res[2]){
      if(GENES.df[GENES.df$GENE==gene,gsub(x = names(which.max(counts)),'*','.',fixed = T)==names(GENES.df)]=='Unk'){
        GENES.df[GENES.df$GENE==gene,gsub(x = names(which.max(counts)),'*','.',fixed = T)==names(GENES.df)] <- strsplit(allele,'*',fixed=T)[[1]][2]
      } else {
        GENES.df[GENES.df$GENE==gene,gsub(x = names(which.max(counts)),'*','.',fixed = T)==names(GENES.df)] <-
          paste(c(GENES.df[GENES.df$GENE==gene,gsub(x = names(which.max(counts)),'*','.',fixed = T)==names(GENES.df)],
                  strsplit(allele,'*',fixed=T)[[1]][2]),collapse=',')
      }
    } else {
      if(GENES.df[GENES.df$GENE==gene,2]=='Unk'){
        GENES.df[GENES.df$GENE==gene,2] <- strsplit(allele,'*',fixed=T)[[1]][2]
      } else {
        GENES.df[GENES.df$GENE==gene,2] <-
          paste(c(GENES.df[GENES.df$GENE==gene,2],
                  strsplit(allele,'*',fixed=T)[[1]][2]),collapse=',')
      }

      if(GENES.df[GENES.df$GENE==gene,3]=='Unk'){
        GENES.df[GENES.df$GENE==gene,3] <- strsplit(allele,'*',fixed=T)[[1]][2]
      } else {
        GENES.df[GENES.df$GENE==gene,3] <-
          paste(c(GENES.df[GENES.df$GENE==gene,3],
                  strsplit(allele,'*',fixed=T)[[1]][2]),collapse=',')
      }

    }



    counts.list[[i]] <- counts
    res.list[[i]] <- res
  }

  len.counts.list <- length(counts.list)
  GENES.df <- cbind(GENES.df,data.frame(ALLELES =  paste(sapply(strsplit(tohap,'*',fixed = T),'[',2),collapse = ','),
                                        PRIORS_ROW = paste(format(HapByPriors,digits = 2),collapse = ','),
                                        PRIORS_COL = paste(format(toHapPriors,digits = 2),collapse = ','),
                                        COUNTS1=paste(counts.list[[1]][order(names(counts.list[[1]])[1:2])],collapse = ','),
                                        MP1=max(res.list[[1]][1:2]),
                                        K1=max(res.list[[1]][1:2])-min(res.list[[1]][1:2]),
                                        ND1=res.list[[1]][3],
                                        COUNTS2=ifelse(length(counts.list) > 1,paste(counts.list[[2]][order(names(counts.list[[2]])[1:2])],collapse = ','),NA),
                                        MP2=ifelse(length(counts.list)>1,max(res.list[[2]][1:2]),NA),
                                        K2=ifelse(length(counts.list)>1,max(res.list[[2]][1:2])-min(res.list[[2]][1:2]),NA),
                                        ND2=ifelse(length(counts.list)>1,res.list[[2]][3],NA),
                                        COUNTS3=ifelse(length(counts.list)>2,paste(counts.list[[3]][order(names(counts.list[[3]])[1:2])],collapse = ','),NA),
                                        MP3=ifelse(length(counts.list)>2,max(res.list[[3]][1:2]),NA),
                                        K3=ifelse(length(counts.list)>2,max(res.list[[3]][1:2])-min(res.list[[3]][1:2]),NA),
                                        ND3=ifelse(length(counts.list)>2,res.list[[3]][3],NA),
                                        COUNTS4=ifelse(length(counts.list)>3,paste(counts.list[[4]][order(names(counts.list[[4]])[1:2])],collapse = ','),NA),
                                        MP4=ifelse(length(counts.list)>3,max(res.list[[4]][1:2]),NA),
                                        K4=ifelse(length(counts.list)>3,max(res.list[[4]][1:2])-min(res.list[[4]][1:2]),NA),
                                        ND4=ifelse(length(counts.list)>3,res.list[[4]][3],NA)
  ))

  return(GENES.df)
}

########################################################################################################
# Haplotype table to plot tables
#
# \code{parseHapTab} Parse the haplotype table for plotting
#
# @param    hap_table             haplotype summary table
# @param    hapBy_alleles         Alleles columns haplotyped by
#
# @return   tables for the three panles in plot
#
parseHapTab <- function(hap_table,hapBy_alleles){

  hap_table <- data.frame(lapply(hap_table, as.character), stringsAsFactors=FALSE)

  # Add count panels
  count.df <- c()
  for(panel in 1:2){
    panel.alleles <- hap_table[[gsub('*','_',hapBy_alleles[panel],fixed = T)]]
    for(i in 1:length(panel.alleles)){
      if(panel.alleles[i]=='Unk'){
        count.df <- rbind(count.df,
                          data.frame(GENE=hap_table$GENE[i],HapBy=hapBy_alleles[panel],
                                     COUNT=0))
      }
      else{
        if(panel.alleles[i]=='Del'){
          count.df <- rbind(count.df,
                            data.frame(GENE=hap_table$GENE[i],HapBy=hapBy_alleles[panel],
                                       COUNT=as.numeric(strsplit(hap_table$COUNTS1[i],',')[[1]][panel])))
        }
        else{
          alleles <- strsplit(panel.alleles[i],',')[[1]]
          for(j in 1:length(alleles)){
            count.df <- rbind(count.df,
                              data.frame(GENE=paste0(hap_table$GENE[i],'*',alleles[j]),HapBy=hapBy_alleles[panel],
                                         COUNT=as.numeric(strsplit(hap_table[i,paste0('COUNTS',j)],',')[[1]][panel])))
          }
        }
      }
    }
  }

  count.df$ALLELES <- sapply(strsplit(as.character(count.df$GENE),'*',fixed = T),'[',2)
  count.df$ALLELES[is.na(count.df$ALLELES)] <- '01' # Mock allele
  count.df$ALLELES <- factor(count.df$ALLELES,levels=c(sort(unique(count.df$ALLELES)),'NA'))
  count.df$GENE <- sapply(strsplit(as.character(count.df$GENE),'*',fixed = T),'[',1)



  # Add K panels
  panel1.alleles <- hap_table[[gsub('*','_',hapBy_alleles[1],fixed = T)]]
  # minimum of Ks if there is more than one allele

  hap_table[is.na(hap_table)] <- 0

  panel1 <- sapply(1:length(panel1.alleles),function(i){
    if(panel1.alleles[i]=='Unk' | panel1.alleles[i]=='Del'){
      min(as.numeric(hap_table[i,paste0('K',1:4)]),na.rm = T)}
    else {min(as.numeric(hap_table[i,paste0('K',match(unlist(strsplit(panel1.alleles[i],',')),
                                                      unlist(strsplit(as.character(hap_table$ALLELES[i]),','))))]))}})
  panel1[panel1==Inf] <- 0

  panel2.alleles <- hap_table[[gsub('*','_',hapBy_alleles[2],fixed = T)]]
  # minimum of Ks if there is more than one allele
  panel2 <- sapply(1:length(panel2.alleles),function(i){
    if(panel2.alleles[i]=='Unk' | panel2.alleles[i]=='Del'){
      min(as.numeric(hap_table[i,paste0('K',1:4)]),na.rm = T)}
    else {min(as.numeric(hap_table[i,paste0('K',match(unlist(strsplit(panel2.alleles[i],',')),
                                                      unlist(strsplit(as.character(hap_table$ALLELES[i]),','))))]))}})
  panel2[panel2==Inf] <- 0

  K.plot <- data.frame(GENE=c(hap_table$GENE,hap_table$GENE),K=c(panel1,panel2),
                       hapBy=c(rep(gsub('IG[H|L|K]',' ',hapBy_alleles[1]),length(panel1)),rep(gsub('IG[H|L|K]',' ',hapBy_alleles[2]),length(panel2))))

  # Add ALLELES
  # Facet by subject
  a <- hap_table[,c('GENE',paste(gsub('*','_',hapBy_alleles[1],fixed = T)))]
  b <- hap_table[,c('GENE',paste(gsub('*','_',hapBy_alleles[2],fixed = T)))]



  a$hapBy <- rep(gsub('*','_',hapBy_alleles[1],fixed = T),nrow(a))
  b$hapBy <- rep(gsub('*','_',hapBy_alleles[2],fixed = T),nrow(b))

  names(a)[2] <- 'ALLELES'
  names(b)[2] <- 'ALLELES'

  geno_sub <-  rbind(a,b)
  geno_sub$GENE <- as.character(geno_sub$GENE)

  names(geno_sub)[1] <- 'GENE'

  return(list(genotype=geno_sub,kval.df=K.plot,count.df=count.df))
}

########################################################################################################
# Binom test for deletion infrence
#
# \code{binom_test_deletion} Infer deletion from binomial test
#
# @param    GENE.usage.df          a data frame of relative gene usage
# @param    cutoff                 a data frame of relative gene usage
# @param    p.val.cutoff           a p value cutoff to detect deletion
# @param    chain                  the IG chain: IGH,IGK,IGL. Default is IGH.
# @param    GENE.loc.IG            the genes by location
#
# @return   data frame with the binomial test results
#

binom_test_deletion <- function(GENE.usage.df,cutoff=0.001,p.val.cutoff=0.01,chain='IGH',GENE.loc.IG)
{
  GENE.usage.df$pval <-  sapply(1:nrow(GENE.usage.df),function(i){if((GENE.usage.df$FRAC[i] < cutoff) &
                                                                                  GENE.usage.df$min_frac[i]!=Inf){
    return(binom.test(x = round(GENE.usage.df$FRAC[i]*GENE.usage.df$NREADS[i]),
                      n = GENE.usage.df$NREADS[i],
                      p = GENE.usage.df$min_frac[i])$p.value)}
    if(GENE.usage.df$min_frac[i]==Inf){
      return(0)
    } else {
      return(1)
    }})

  ### P.binom to detect deletion or cnv
  GENE.usage.df$foradj <-  sapply(1:nrow(GENE.usage.df),function(i){
    if(GENE.usage.df$FRAC[i] < cutoff & GENE.usage.df$min_frac[i]!=Inf){return(paste0(GENE.usage.df$GENE[i],'_',0))};
    if(GENE.usage.df$min_frac[i]==Inf) {
      return(paste0(GENE.usage.df$GENE[i],'_',1))
    } else {
      return(paste0(GENE.usage.df$GENE[i],'_',2))
    }
  })
  GENE.usage.df <- GENE.usage.df %>% group_by(foradj) %>% mutate(pval_adj=p.adjust(pval, method = 'BH'))



  GENE.usage.df$col <- sapply(1:nrow(GENE.usage.df),function(i){if(GENE.usage.df$pval_adj[i] <= p.val.cutoff){
    if(GENE.usage.df$FRAC[i] < cutoff & GENE.usage.df$min_frac[i]!=Inf){return('Deletion')};
    if(GENE.usage.df$min_frac[i]==Inf) {
      return('NA')
    }} else {
      return('No Deletion')
    }
  })


  GENE.usage.df$GENE <- factor(GENE.usage.df$GENE, levels=GENE.loc.IG)
  GENE.usage.df$col <- factor(GENE.usage.df$col, levels=c('Deletion','No Deletion','NA') )

  return(GENE.usage.df)
}
