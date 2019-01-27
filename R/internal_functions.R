# Internal functions -----------------------------------------------------

#' @include rabhit.R
NULL

########################################################################################################
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
  GENES.df.num <- reshape2::melt(df)

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
# \code{parseHapTab} Parse the haplotype table for each panel in the haplotype plot
#
# @param    hap_table             haplotype summary table
# @param    chain                 the Ig chain: IGH,IGK,IGL. Default is IGH.
# @param    hapBy_alleles         Alleles columns haplotyped by
#
# @return   list of data frames for plotting
#
parseHapTab <- function(hap_table, chain = c('IGH','IGK','IGL'),
                        df_ToReturn = c('all','geno.df','kval.df','count.df')){

  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(missing(df_ToReturn)) {
    df_ToReturn='all'
  }

  hap_table <- data.frame(lapply(hap_table, as.character), stringsAsFactors=FALSE)
  hapBy_cols = names(hap_table)[grep(chain,names(hap_table))]
  hapBy_alleles = gsub('_','*',hapBy_cols)



  # Count panel data frame
  count.df <- c()
  for(panel in 1:2){
    panel.alleles <- hap_table[[hapBy_cols[panel]]]
    for(i in 1:length(panel.alleles)){
      if(panel.alleles[i]=='Unk' | panel.alleles[i]=='NR'){
        count.df <- rbind(count.df,
                          data.frame(GENE=hap_table$GENE[i],hapBy=hapBy_alleles[panel],
                                     COUNT=0))
      }
      else{
        if(panel.alleles[i]=='Del'){
          count.df <- rbind(count.df,
                            data.frame(GENE=hap_table$GENE[i],hapBy=hapBy_alleles[panel],
                                       COUNT=as.numeric(strsplit(hap_table$COUNTS1[i],',')[[1]][panel])))
        }
        else{
          alleles <- strsplit(panel.alleles[i],',')[[1]]
          for(j in 1:length(alleles)){
            count.df <- rbind(count.df,
                              data.frame(GENE=paste0(hap_table$GENE[i],'*',alleles[j]),
                                         hapBy=hapBy_alleles[panel],
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

  ## TO visualy make coutns of 1 not look like 0 , one is added
  count.df$COUNT2 <- ifelse(count.df$hapBy == hapBy_alleles[1],-1*log10(as.numeric(count.df$COUNT)+1),log10(as.numeric(count.df$COUNT)+1))
  count.df$COUNT2[count.df$COUNT2==Inf |count.df$COUNT2== -Inf  ] <- 0


  # K values panels data frame
  panel1.alleles <- hap_table[[hapBy_cols[1]]]
  # minimum of Ks if there is more than one allele

  hap_table[is.na(hap_table)] <- 0

  panel1 <- sapply(1:length(panel1.alleles),function(i){
    if(panel1.alleles[i]=='Unk' | panel1.alleles[i]=='Del' | panel1.alleles[i]=='NR'){
      min(as.numeric(hap_table[i,paste0('K',1:4)]),na.rm = T)}
    else {min(as.numeric(hap_table[i,paste0('K',match(unlist(strsplit(panel1.alleles[i],',')),
                                                      unlist(strsplit(as.character(hap_table$ALLELES[i]),','))))]),na.rm = T)}})
  panel1[panel1==Inf] <- 0

  panel2.alleles <- hap_table[[hapBy_cols[2]]]
  # minimum of Ks if there is more than one allele
  panel2 <- sapply(1:length(panel2.alleles),function(i){
    if(panel2.alleles[i]=='Unk' | panel2.alleles[i]=='Del' | panel2.alleles[i]=='NR'){
      min(as.numeric(hap_table[i,paste0('K',1:4)]),na.rm = T)}
    else {min(as.numeric(hap_table[i,paste0('K',match(unlist(strsplit(panel2.alleles[i],',')),
                                                      unlist(strsplit(as.character(hap_table$ALLELES[i]),','))))]))}})
  panel2[panel2==Inf] <- 0

  kval.df <- data.frame(GENE=c(hap_table$GENE,hap_table$GENE),K=c(panel1,panel2),
                        hapBy=c(rep(hapBy_alleles[1],length(panel1)),
                                rep(hapBy_alleles[2],length(panel2))))

  kval.df$K_GROUPED <- bin_data(kval.df$K, bins=c(0, 1,2,3,4,5,10,20,50,Inf), binType = "explicit")

  # Alleles panel data frame
  a <- hap_table[,c('GENE',hapBy_cols[1])]
  b <- hap_table[,c('GENE',hapBy_cols[2])]
  a$hapBy <- rep(hapBy_cols[1],nrow(a))
  b$hapBy <- rep(hapBy_cols[2],nrow(b))
  names(a)[2] <- 'ALLELES'
  names(b)[2] <- 'ALLELES'
  genotype <-  rbind(a,b)

  alleles = strsplit(genotype$ALLELES, ",")
  geno.df = genotype
  r = 1
  for (g in 1:nrow(genotype)) {
    for (a in 1:length(alleles[[g]])) {
      geno.df[r, ] = genotype[g, ]
      geno.df[r, ]$ALLELES = alleles[[g]][a]
      r = r + 1
    }
  }
  geno.df$hapBy <- gsub('_','*',geno.df$hapBy,fixed = T)

  parsed_hap_table <- list(geno.df=geno.df,kval.df=kval.df,count.df=count.df)
  if(df_ToReturn=='all') return(parsed_hap_table)
  else parsed_hap_table[[df_ToReturn]]
}

########################################################################################################
# Sort data frame by genes
#
# \code{sortDFByGene} Sort the \code{data.frame} by the genes names or position. For sorting by gene names the \code{sortAlleles} function by TIgGER is used.
# For sorting by position the defualt package gene location list is used.
#
# @param    DATA                 data frame to sort
# @param    chain                the Ig chain: IGH,IGK,IGL. Default is IGH.
# @param    method               the method for sorting the genes. If by 'name' the genes in the output are ordered lexicographically,
# if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
# @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#
# @return   sorted \code{data.frame}
#
sortDFByGene <- function(DATA, chain = c('IGH','IGK','IGL'), method = c("name", "position"), removeIGH = FALSE){
  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(missing(method)) {
    method='position'
  }
  method <- match.arg(method)

  if(method=='name'){
    DATA$GENE <- factor(DATA$GENE, levels = rev(sortAlleles(unique(DATA$GENE), method = method)))
    if(removeIGH){
      DATA$GENE <- gsub('IG[H|K|L]','',DATA$GENE)
      DATA$GENE <- factor(DATA$GENE, levels = rev(sortAlleles(unique(DATA$GENE), method = method)))
      DATA$hapBy <- gsub('IG[H|K|L]','',DATA$hapBy)
    }
  } else {
    GENE.loc.tmp <- GENE.loc[[chain]]
    names(GENE.loc.tmp) <- GENE.loc.tmp
    DATA$GENE <- factor(DATA$GENE, levels = rev(GENE.loc.tmp))
    if(removeIGH){
      GENE.loc.tmp <- gsub('IG[H|K|L]','',GENE.loc.tmp)
      names(GENE.loc.tmp) <- GENE.loc.tmp
      DATA$GENE <- gsub('IG[H|K|L]','',DATA$GENE)
      DATA$GENE <- factor(DATA$GENE, levels = rev(GENE.loc.tmp))
      DATA$hapBy <- gsub('IG[H|K|L]','',DATA$hapBy)
    }
  }

  return(DATA)
}

########################################################################################################
# Calculate Jaacardian distance for haplotypes
#
# \code{calcJacc} Takes as an input two haplotypes and calculates the Jaacardian distance.
#
# @param   vec1A           chromosome A haplotype for first individual.
# @param   vec1B           chromosome B haplotype for first individual.
# @param   vec2A           chromosome A haplotype for second individual.
# @param   vec2B           chromosome B haplotype for second individual.
# @param   method          the method to be used for calculating. pooled - All alleles and all genes taken together (assuming that all genes appear and ordered the same in both vectors)
# geneByGene - For each gene separetly and then calculates average distance.
# @param  naRm             if 'TRUE' ingnores genes from both samples in which there is no data, else the Jaccardian distance for those gene defined as 1 (i.e the maximal distance)
# If both are NAs will remove it either way.
# @param   Kweight         if 'TRUE' the Jaacardian distance is weighted by the K value of the haplotype.
# @param   k1A             chromosome A haplotype K value for first individual.
# @param   k1B             chromosome B haplotype K value for first individual.
# @param   k2A             chromosome A haplotype K value for second individual.
# @param   k2B             chromosome B haplotype K value for second individual.
#
# @return   Jaacardian distance value
#
distJACC <- function(vecA, vecB, naRm = TRUE){
  if((!is.na(vecA) & !is.na(vecB)) & (vecA != 'Unk' & vecB != 'Unk')){
    v1 <- unlist(strsplit(vecA,','))
    v2 <- unlist(strsplit(vecB,','))

    if((any(grepl('_[0-9][0-9]',v1))|any(grepl('_[0-9][0-9]',v2)))&(length(intersect(v1,v2))!=length(unique(c(v1,v2))))){
      inter = 0
      tot = length((c(v1,v2)))
      for(i in v1){
        v1_n <- unlist(strsplit(i,'_'))
        for(ii in v2){
          v2_n <- unlist(strsplit(ii,'_'))
          if(length(intersect(v1_n,v2_n))!=0){
            inter = inter + 1
            tot = tot - 1
          }

        }

      }
      dist <- inter/tot
    }
    else dist <- length(intersect(v1,v2))/length(unique(c(v1,v2)))
  } else {
    if(naRm){
      dist <- NA
    } else {
      if((is.na(vecA) & is.na(vecB)) || (vecA == 'Unk' & vecB == 'Unk')){
        dist <- NA
      }
      else {
        dist <- 0

      }
    }

  }

}

calcJacc <- function(vec1A, vec1B, vec2A, vec2B, method=c('pooled', 'geneByGene'), naRm = TRUE, Kweight = FALSE, k1A, k1B, k2A, k2B, special_jacc = FALSE){
  vec1A <- as.character(vec1A)
  vec1B <- as.character(vec1B)
  vec2A <- as.character(vec2A)
  vec2B <- as.character(vec2B)
  if(method=='geneByGene'){
    # Calculate Jaccardian distance for each gene and then
    jacc <- c()
    jacc <- sapply(1:length(vec1A),function(i){
      distA <- distJACC(vec1A[i],vec2A[i], naRm)
      distB <- distJACC(vec1B[i],vec2B[i], naRm)
      jacc <- c(jacc,rowMeans(data.frame(distA, distB), na.rm = TRUE))
    })

    if(Kweight){
      k1A[k1A==Inf] <- 0;k2A[k2A==Inf] <- 0;k1B[k1B==Inf] <- 0;k2B[k2B==Inf] <- 0
      KavgA <- sapply(1:length(k1A),function(x){mean(c(k1A[x],k2A[x]),na.rm = T)})
      if(is.list(KavgA)){KavgA[sapply(KavgA, is.null)] <- NA; KavgA <- unlist(KavgA)}
      KavgB <- sapply(1:length(k1B),function(x){mean(c(k1B[x],k2B[x]),na.rm = T)})
      if(is.list(KavgB)){KavgB[sapply(KavgB, is.null)] <- NA; KavgB <- unlist(KavgB)}

      Kavg <- rowMeans(cbind(KavgA,KavgB),na.rm = T)
      jacc <- weighted.mean(jacc,Kavg,na.rm = T)
      return(1-jacc)
    }
    return(mean(1-jacc,na.rm = T))
  }
  if(method=='pooled'){

    v1 <- unlist(sapply(1:length(vec1A), function(x){paste(paste0('G',x),unlist(strsplit(vec1A[[x]],',')),sep = '_')}))
    v1 <- c(v1,unlist(sapply(1:length(vec1B), function(x){paste(paste0('G',x),unlist(strsplit(vec1B[[x]],',')),sep = '_')})))
    v2 <- unlist(sapply(1:length(vec2A), function(x){paste(paste0('G',x),unlist(strsplit(vec2A[[x]],',')),sep = '_')}))
    v2 <- c(v2,unlist(sapply(1:length(vec2B), function(x){paste(paste0('G',x),unlist(strsplit(vec2B[[x]],',')),sep = '_')})))

    v1 <- v1[-grep('NA',v1,fixed=T)]
    v2 <- v2[-grep('NA',v2,fixed=T)]

    jacc <- length(intersect(v1,v2))/length(unique(c(v1,v2)))
    return(1-jacc)
  }
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

########################################################################################################
# Creates the allele color palette for haplotype graphical output
#
# \code{alleleHapPalette} Takes a list of the haplotype alleles and returns the allele color palette.
#
# @param   hap_alleles          a list of the haplotype alleles.
#
# @return   Haplotype allele color palette
#
alleleHapPalette <- function(hap_alleles){

  AlleleCol <- grep('[012]',unique(hap_alleles),value = T,perl = T)
  AlleleCol.tmp <- sort(unique(sapply(strsplit(AlleleCol,'_'),'[',1)))
  tmp.col <- ALLELE_PALETTE[AlleleCol.tmp]

  novels <- grep('_',AlleleCol,value = T)
  if(length(novels) > 0){
    novels.col <- ALLELE_PALETTE[sapply(strsplit(novels,'_'),'[',1)]
    names(novels.col) <- novels
    alleles.comb <- c(tmp.col,novels.col)[order(names(c(tmp.col,novels.col)))]
  } else {
    alleles.comb <- c(tmp.col)[order(names(c(tmp.col)))]

  }

  AlleleCol<- names(c(alleles.comb,Unk='#dedede',Del='#6d6d6d',NR='#000000'))
  names(AlleleCol) <- c(alleles.comb,Unk='#dedede',Del='#6d6d6d',NR='#000000')

  transper <- sapply(AlleleCol,function(x){
    if(grepl('_',x)){
      mom_allele <- strsplit(x,'_')[[1]][1];
      all_novel <- grep(paste0(mom_allele,'_'),AlleleCol,value=T);
      if(length(all_novel)==1){return(0.5)};
      if(length(all_novel)==2){m=which(all_novel==x);return(ifelse(m==1,0.6,0.3))}
      if(length(all_novel)==3){m=which(all_novel==x);if(m==1){return(0.6)} ; return(ifelse(m==2,0.4,0.2))}
      if(length(all_novel)>9){m=which(all_novel==x);if(m==1){return(1)} ; return(1-m/20)}
      if(length(all_novel)>3){m=which(all_novel==x);if(m==1){return(0.85)} ; return(0.85-m/10)}
    }else (1)})
  names(transper) <- AlleleCol

  #remove 'mother' allele if added (when there is no germline allele but there is a novel)
  AlleleCol <- AlleleCol[AlleleCol %in% c(sort(grep('[012]',unique(hap_alleles),value = T,perl = T)),'Unk','Del','NR')]

  transper <- transper[names(transper) %in% AlleleCol ]

  return(list(transper=transper,AlleleCol=AlleleCol))
}

########################################################################################################
#' Detect non reliable gene assignment
#'
#' \code{nonReliableGenes} Takes a \code{data.frame} in Change-O format  and detect non reliable IGHV genes. A non reliable gene is
#' when the ratio of the multiple assignment with a gene is below the threshold.
#'
#'  @details Function accepts a \code{data.frame} in Change-O format (\link{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{"SUBJECT"}: subject names
#'   \item \code{"V_CALL"}: V allele call(s) (in an IMGT format)
#'   \item \code{"D_CALL"}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{"J_CALL"}: J allele call(s) (in an IMGT format)
#' }
#'
#' @param   clip_db            a \code{data.frame} in Change-O format. See details.
#' @param   thresh             the threshold to consider non reliable gene. Defualt is 0.9
#' @param   chain              the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param   appearance         the minimun fraction of gene appearance to be considered for reliability check. Defualt is 0.01.
#' @return  a nested list of non reliable genes for all subject.
#'
#' @examples
#' # Example IGHV call list
#' V_CALL <- c('IGHV1-69*01','IGHV1-69*01','IGHV1-69*01,IGHV1-69*02','IGHV4-59*01,IGHV4-61*01','IGHV4-59*01,IGHV4-31*02','IGHV4-59*01')
#'
#' # Detect non reliable genes
#' nonReliableGenes(V_CALL)
#'
#' @export

nonReliableGenes <- function(clip_db,thresh=0.9,appearance=0.01,chain=c('IGH','IGK','IGL')){

  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(!("SUBJECT" %in% names(clip_db))){clip_db$SUBJECT <- 'S1'}
  GENE.loc.NoPseudo <- GENE.loc[[chain]][!grepl('OR',GENE.loc[[chain]])]
  GENE.loc.NoPseudo <- GENE.loc.NoPseudo[!grepl('NL',GENE.loc.NoPseudo)]
  GENE.loc.NoPseudo <- GENE.loc.NoPseudo[!(GENE.loc.NoPseudo %in% PSEUDO[[chain]])]
  GENE.loc.NoPseudo <- GENE.loc.NoPseudo[grepl('V',GENE.loc.NoPseudo)]

  non_reliable_genes <- c()
  for(sample_name in unique(clip_db$SUBJECT)){
    clip_db_sub = clip_db[clip_db$SUBJECT==sample_name,]
    gene_call <- getGene(clip_db_sub$V_CALL,first = F,strip_d = F)
    NR_tmp <- c()
    for(gene in GENE.loc.NoPseudo){
      sa <- length(grep(paste0('^',gene,'$'),gene_call))
      total <- length(grep(paste0('^',gene,',|',',',gene,'$|',',',gene,',|','^',gene,'$'),gene_call))

      if(total/length(gene_call)>=appearance) if((sa/total)<thresh) NR_tmp <- c(NR_tmp,gene)
    }

    if(!is.null(NR_tmp)) non_reliable_genes[[sample_name]] <- NR_tmp
  }
  return(non_reliable_genes)
}
