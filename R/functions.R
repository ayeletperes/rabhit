# RAbHIT functions -----------------------------------------------------

#' @include rabhit.R
#' @include internal_functions.R
NULL

#' Anchor gene haplotype inference
#'
#' \code{createFullHaplotype} infers haplotype based on an anchor gene.
#'
#' @details Function accepts a \code{data.frame} in Change-O format (\link{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{"SUBJECT"}: subject names
#'   \item \code{"V_CALL"}: V allele call(s) (in an IMGT format)
#'   \item \code{"D_CALL"}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{"J_CALL"}: J allele call(s) (in an IMGT format)
#' }
#'
#' @param    clip_db               a \code{data.frame} in Change-O format. See details.
#' @param    toHap_col             a vector of column names for which a haplotype should be inferred. Default is V_CALL and D_CALL.
#' @param    hapBy_col             column name of the anchor gene. Default is J_CALL.
#' @param    hapBy                 a string of the anchor gene name. Default is IGHJ6.
#' @param    toHap_GERM            a vector of named nucleotide germline sequences matching the allele calls in \code{toHap_col} columns in clip_db.
#' @param    relative_freq_priors  if TRUE, the priors for Bayesian inference are estimated from the relative frequencies in clip_db. Else, priors are set to \code{c(0.5,0.5)}. Defualt is TRUE
#' @param    kThreshDel            The minimum lK (log10 of the Bayes factor) to call a deletion. Defualt is 3.
#' @param    rmPseudo              if TRUE non-functional and pseudo genes are removed. Defualt is TRUE.
#' @param    deleted_genes         Double chromosome deletion summary table. A \code{data.frame} created by \code{binom_test_deletion}.
#' @param    min_minor_fraction    the minimum minor allele fraction to be used as an anchor gene. Default is 0.3
#' @param    chain                 the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    supress_print         if TRUE the function does not print summary. Dufault FALSE
#'
#' @return   a haplotype summary table
#'
#'
#' @examples
#' # Load example data and germlines
#' data(sample_db, HVGERM, HDGERM)
#'
#' # Infering haplotype
#' hap_df = createFullHaplotype(sample_db,toHap_col=c("V_CALL","D_CALL"),hapBy_col="J_CALL",hapBy="IGHJ6",toHap_GERM=c(HVGERM,HDGERM));head(hap_df);
#'
#' @export

createFullHaplotype <- function(clip_db,toHap_col=c("V_CALL","D_CALL"),hapBy_col="J_CALL",hapBy="IGHJ6",
                              toHap_GERM,relative_freq_priors=TRUE,kThreshDel=3,rmPseudo=TRUE,deleted_genes=c(),
                              min_minor_fraction=0.3,chain=c('IGH','IGK','IGL'),supress_print=FALSE){
  if(missing(chain)) {
        chain='IGH'
    }
  chain <- match.arg(chain)

  if(!("SUBJECT" %in% names(clip_db))){clip_db$SUBJECT <- rep('S1',nrow(clip_db))}


  PSEUDO <- PSEUDO[[chain]]
  GENE.loc <- GENE.loc[[chain]]

  haplo_db <- c()
  for(sample_name in unique(clip_db$SUBJECT)){
      clip_db_sub = clip_db[clip_db$SUBJECT==sample_name,]
      if(is.data.frame(deleted_genes)) deleted_genes_vec <- deleted_genes %>% filter(SUBJECT==sample_name,DELETION=='Deletion') %>% select(GENE) %>% pull()
      else deleted_genes_vec <-c()
      #Number of iniial sequences
      nrows1 <- nrow(clip_db_sub)
      ### Check if haplotype can be infered by the specific gene in the data set.
      ### Only relevant genes with one assignment
      IND <- apply(clip_db_sub,1,function(x){sum(grep(',',x[c(toHap_col,hapBy_col)],invert = F))})
      IND <- names(IND[IND==0])
      clip_db_sub <- clip_db_sub[IND,]

      #Number of post multiple assignment reduction sequences
      nrows2 <- nrow(clip_db_sub)

      if(!supress_print) cat(paste0('In sample ',sample_name,', ',nrows1-nrows2,' sequnces were removed due to multiple assignments,\n ',nrows2,' sequences left.\n'))

      hapBy_alleles <- sort(unique(grep(pattern = paste0(hapBy,'*'),x = clip_db_sub[,paste(hapBy_col)],value = T,fixed = T)))

      hapBy_alleles_table <- table(grep(pattern = paste0(hapBy,'*'),x = clip_db_sub[,paste(hapBy_col)],value = T,fixed = T))
      if(length(unique(hapBy_alleles))!=2) stop("Can not haplotype by more or less than two alleles")
      if(min(hapBy_alleles_table)/sum(hapBy_alleles_table) < min_minor_fraction) stop("Can not haplotype, minor allele fraction lower than the cutoff set by the user")


      for(i in 1:length(toHap_col)){
        ## Remove "none" in D assignments
        if(i > 1){
          GENES <- c(GENES,unique(sapply(strsplit(clip_db_sub[,paste(toHap_col[i])][clip_db_sub[,paste(hapBy_col)] %in% hapBy_alleles & grepl('IG',clip_db_sub[,paste(toHap_col[i])]) ],'*',fixed = T),'[',1)))
        } else {
          GENES <- unique(sapply(strsplit(clip_db_sub[,paste(toHap_col[i])][clip_db_sub[,paste(hapBy_col)] %in% hapBy_alleles & grepl('IG',clip_db_sub[,paste(toHap_col[i])]) ],'*',fixed = T),'[',1))

        }

      }


      GENES.ref <- unique(sapply(strsplit(names(toHap_GERM),'*',fixed = T),'[',1))

      GENES.df <- data.frame(GENE=GENES.ref,g1=rep('Unk',length(GENES.ref)),g2=rep('Unk',length(GENES.ref)),stringsAsFactors = F)
      names(GENES.df)[2:3] <- sort(hapBy_alleles)


      if(rmPseudo){
        GENES <- GENES[!grepl('OR',GENES)]
        GENES <- GENES[!grepl('NL',GENES)]
        GENES <- GENES[!(GENES %in% PSEUDO)]

        GENES.ref <- GENES.ref[!grepl('OR',GENES.ref)]
        GENES.ref <- GENES.ref[!grepl('NL',GENES.ref)]
        GENES.ref <- GENES.ref[!(GENES.ref %in% PSEUDO)]

        GENES.df <- GENES.df[!grepl('OR',GENES.df$GENE),]
        GENES.df <- GENES.df[!grepl('NL',GENES.df$GENE),]
        GENES.df <- GENES.df[!(GENES.df$GENE %in% PSEUDO),]
      }

      GENES.df.num <- c()
      GENES.df.num_ToPlot <- c()


      for(G in intersect(GENES,GENES.ref)){


        if(grepl(paste0(chain,'V'),G)) toHap_col_tmp <- 'V_CALL'
        if(grepl(paste0(chain,'D'),G)) toHap_col_tmp <- 'D_CALL'
        if(grepl(paste0(chain,'J'),G)) toHap_col_tmp <- 'J_CALL'

        clip_db_sub.G <- clip_db_sub[grepl(paste0(G,"*"),clip_db_sub[,toHap_col_tmp],fixed = T) & (clip_db_sub[,hapBy_col] %in% hapBy_alleles),]

        tmp <- table(clip_db_sub.G[,toHap_col_tmp],clip_db_sub.G[,hapBy_col])
        relFreq <- min(rowSums(tmp))/sum(rowSums(tmp))
        GENES.df.num_ToPlot <- rbind(GENES.df.num_ToPlot,reshape2::melt(tmp))
        # if one column add the second
        if(ncol(tmp)==1){
          toadd <- setdiff(hapBy_alleles,colnames(tmp))
          tmp <- cbind(tmp,rep(0,nrow(tmp)))
          colnames(tmp)[2] <- toadd

          coltmp <- colnames(tmp)
          rowtmp <- rownames(tmp)

          tmp <- as.data.frame(tmp)
          tmp <- tmp[order(colnames(tmp))]
          tmp <- as.matrix(tmp)
        }


        if(relative_freq_priors){

          clip_db_sub.toHap <- clip_db_sub[clip_db_sub[,hapBy_col] %in% hapBy_alleles,]
          hapBy_priors <- table(clip_db_sub.toHap[,hapBy_col])/sum(table(clip_db_sub.toHap[,hapBy_col]))

          clip_db_sub.hapBy <- clip_db_sub[(grepl(paste0(G,"*"),clip_db_sub[,toHap_col_tmp],fixed = T)) ,]
          clip_db_sub.hapBy <- clip_db_sub.hapBy[clip_db_sub.hapBy[,toHap_col_tmp] %in% rownames(tmp), ]
          toHap_priors <- table(clip_db_sub.hapBy[,toHap_col_tmp])/sum(table(clip_db_sub.hapBy[,toHap_col_tmp]))

          if(length(toHap_priors) != nrow(tmp) ){
            toHap_priors_tmp <- c(rep(0,nrow(tmp)))
            names(toHap_priors_tmp) <- rownames(tmp)
            for(i in names(toHap_priors)){
              toHap_priors_tmp[i] <- toHap_priors[i]
            }

            toHap_priors <- toHap_priors_tmp
          }

          hap.df <- createHaplotypeTable(tmp,HapByPriors = hapBy_priors,toHapByCol = TRUE ,toHapPriors = toHap_priors)

          relFreqDf.tmp <- data.frame(c(sample_name,G,format(relFreq,digits = 3),sum(grepl(',',hap.df[1,2:3],fixed = T)),hap.df[,2:length(hap.df)]),stringsAsFactors = F)
          names(relFreqDf.tmp)[1:6]<- c('SUBJECT','GENE','MinorFraction','DoubleAllele',gsub(pattern = '*','_',hapBy_alleles,fixed = T))
          GENES.df.num <- rbind(GENES.df.num,relFreqDf.tmp)
        } else {
          hap.df <- createHaplotypeTable(tmp)
          relFreqDf.tmp <- data.frame(c(sample_name,G,format(relFreq,digits = 3),sum(grepl(',',hap.df[1,2:3],fixed = T)),hap.df[,2:length(hap.df)]),stringsAsFactors = F)
          names(relFreqDf.tmp)[1:6]<- c('SUBJECT','GENE','MinorFraction','DoubleAllele',gsub(pattern = '*','_',hapBy_alleles,fixed = T))
          GENES.df.num <- rbind(GENES.df.num,relFreqDf.tmp)

        }
      }





      ## Fill deleted according to a k thershold
      delIDX <-  which(GENES.df.num[paste(gsub('*','_',hapBy_alleles[1],fixed = T))][,1] == 'Unk' & sapply(1:nrow(GENES.df.num),function(x){min(GENES.df.num[x,grepl('K[1-9]',names(GENES.df.num))],na.rm = T)}) >= kThreshDel)

      GENES.df.num[delIDX,paste(gsub('*','_',hapBy_alleles[1],fixed = T))] <- "Del"

      delIDX <-  which(GENES.df.num[paste(gsub('*','_',hapBy_alleles[2],fixed = T))][,1] == 'Unk' & sapply(1:nrow(GENES.df.num),function(x){min(GENES.df.num[x,grepl('K[1-9]',names(GENES.df.num))],na.rm = T)}) >= kThreshDel)

      GENES.df.num[delIDX,paste(gsub('*','_',hapBy_alleles[2],fixed = T))] <- "Del"

      ## Add as iunknown genes that do not appear in the individual and mark them as unknown

      GENES.MISSING <- GENES.df$GENE[!(GENES.df$GENE %in% GENES.df.num$GENE)]
      if(length(GENES.MISSING) > 0){
        m <- length(GENES.MISSING)

        sub.df <- data.frame(do.call(rbind,lapply(1:m, function(x){c(sample_name,NA,0,0,'Unk','Unk',rep(NA,ncol(GENES.df.num)-6))})))
        names(sub.df) <- names(GENES.df.num)

        sub.df$GENE <- GENES.MISSING


        GENES.df.num <- rbind(GENES.df.num,sub.df)}

      #### Remove deleted genes according to known list

      if(length(deleted_genes_vec)!=0){
        GENES.df.num[,5] <- sapply(1:nrow(GENES.df.num),function(i){if(GENES.df.num$GENE[i] %in% deleted_genes_vec){
          return('Del')
        } else {
          return( GENES.df.num[i,5])
        }})
        GENES.df.num[,6] <- sapply(1:nrow(GENES.df.num),function(i){if(GENES.df.num$GENE[i] %in% deleted_genes_vec){
          return('Del')
        } else {
          return( GENES.df.num[i,6])
        }
        })
        GENES.df.num['K1'] <- sapply(1:nrow(GENES.df.num),function(i){if(GENES.df.num$GENE[i] %in% deleted_genes_vec){
          return(NA)
        } else {
          return( GENES.df.num[i,'K1'])
        }
        })
      }
    haplo_db <- rbind(haplo_db,GENES.df.num)
  }

  return(haplo_db)

}

##########################################################################
#' Double chromosome deletion by relative gene usage
#'
#' \code{deletionsByBinom} detects chromosomal deletion.
#'
#' @details Function accepts a \code{data.frame} in Change-O format (\link{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{"SUBJECT"}: subject names
#'   \item \code{"V_CALL"}: V allele call(s) (in an IMGT format)
#'   \item \code{"D_CALL"}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{"J_CALL"}: J allele call(s) (in an IMGT format)
#' }
#'
#' @param    clip_db               a \code{data.frame} in Change-O format. See details.
#' @param    chain                 the IG chain: IGH,IGK,IGL. Default is IGH.
#' @return  data frame with double chromosome gene deletions
#'
#'
#' @examples
#' # Load example data and germlines
#' data(sample_db)
#'
#' # Infering haplotype
#' del_binom_df = deletionsByBinom(sample_db);head(del_binom_df);
#'
#' @export

deletionsByBinom <- function(clip_db,chain=c('IGH','IGK','IGL')){
  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(!("SUBJECT" %in% names(clip_db))){clip_db$SUBJECT <- rep('S1',nrow(clip_db))}


  PSEUDO <- PSEUDO[[chain]]
  GENE.loc <- GENE.loc[[chain]]

  GENE.loc.NoPseudo <- GENE.loc[!grepl('OR',GENE.loc)]
  GENE.loc.NoPseudo <- GENE.loc.NoPseudo[!grepl('NL',GENE.loc.NoPseudo)]
  GENE.loc.NoPseudo <- GENE.loc.NoPseudo[!(GENE.loc.NoPseudo %in% PSEUDO)]

  GENE.usage <-vector("list", length = 80)
  names(GENE.usage) <- GENE.loc.NoPseudo
  for(samp in unique(clip_db$SUBJECT)){
    clip_db_sub <- clip_db[clip_db$SUBJECT==samp,]
    # V gene distribution
    V_CALLS <- table(sapply(strsplit(clip_db_sub$V_CALL,'*',fixed = T),'[',1))
    V_CALLS_freq <- V_CALLS/sum(V_CALLS)
    for(v in GENE.loc){
      if(v %in% names(V_CALLS)){
        GENE.usage[[v]] <- c( GENE.usage[[v]],V_CALLS_freq[v])
        names(GENE.usage[[v]])[length(GENE.usage[[v]])] <- samp
      } else {
        if(grepl(paste0(chain,'V'),v)){
          GENE.usage[[v]] <- c( GENE.usage[[v]],0)
          names(GENE.usage[[v]])[length(GENE.usage[[v]])] <- samp
        }
      }
    }
    # D gene distribution
    if(chain=='IGH'){
      D_SINGLE <- grep(pattern = ',',clip_db_sub$D_CALL,invert = T)
      D_CALLS <- table(sapply(strsplit(clip_db_sub$D_CALL[D_SINGLE],'*',fixed = T),'[',1))
      D_CALLS_freq <- D_CALLS/sum(D_CALLS)
      for(d in GENE.loc){
        if(d %in% names(D_CALLS)){
          GENE.usage[[d]] <- c( GENE.usage[[d]],D_CALLS_freq[d])
          names(GENE.usage[[d]])[length(GENE.usage[[d]])] <- samp
        } else {
          if(grepl(paste0(chain,'D'),d)){
            GENE.usage[[d]] <- c( GENE.usage[[d]],0)
            names(GENE.usage[[d]])[length(GENE.usage[[d]])] <- samp
          }
        }
      }}
    # J gene distribution
    J_CALLS <- table(sapply(strsplit(clip_db_sub$J_CALL,'*',fixed = T),'[',1))
    J_CALLS_freq <- J_CALLS/sum(J_CALLS)
    for(j in GENE.loc){
      if(j %in% names(J_CALLS)){
        GENE.usage[[j]] <- c( GENE.usage[[j]],J_CALLS_freq[j])
        names(GENE.usage[[j]])[length(GENE.usage[[j]])] <- samp
      } else {
        if(grepl(paste0(chain,'J'),j)){
          GENE.usage[[j]] <- c( GENE.usage[[j]],0)
          names(GENE.usage[[j]])[length(GENE.usage[[j]])] <- samp
        }
      }
    }
  }

  GENE.usage.mean.dist <- sapply(GENE.usage,function(x){(x-mean(x))/sd(x)})
  GENE.usage.mean.dist.mlt <- reshape2::melt(GENE.usage.mean.dist)
  no.existent.genes<- which(colSums(GENE.usage.mean.dist,na.rm = T) == 0)
  GENE.usage.mean.dist.mlt <- GENE.usage.mean.dist.mlt[!(GENE.usage.mean.dist.mlt$Var2 %in% names(no.existent.genes)), ]
  GENE.usage.mean.dist.mlt$value1 <- cut(GENE.usage.mean.dist.mlt$value,breaks = c(-Inf,-3:3,Inf),right = FALSE)


  ###  Gene usage for violin plot for paper
  SAMPLE.SIZE.V <- sapply(unique(clip_db$SUBJECT),function(x){sub <- clip_db[clip_db$SUBJECT==x,];nrow(sub[!grepl(',',sub$V_CALL),])})
  if(chain=='IGH'){SAMPLE.SIZE.D <- sapply(unique(clip_db$SUBJECT),function(x){sub <- clip_db[clip_db$SUBJECT==x,];nrow(sub[!grepl(',',sub$D_CALL),])})}
  SAMPLE.SIZE.J <- sapply(unique(clip_db$SUBJECT),function(x){sub <- clip_db[clip_db$SUBJECT==x,];nrow(sub[!grepl(',',sub$J_CALL),])})
  SAMPLE.SIZE <- sapply(unique(clip_db$SUBJECT),function(x) nrow(clip_db[clip_db$SUBJECT==x,]))

  gusage <- unlist(GENE.usage)
  gusage.gene <- sapply(strsplit(names(unlist(GENE.usage)),'.',fixed = T),'[',1)
  gusage.samp <- sapply(strsplit(names(unlist(GENE.usage)),'.',fixed = T),'[',2)
  GENE.usage.df <- data.frame(SUBJECT=gusage.samp,GENE=gusage.gene,FRAC=gusage,stringsAsFactors = F,row.names = NULL)

  GENE.usage.df <- GENE.usage.df %>% filter(GENE %in% GENE.loc.NoPseudo)
  GENE.usage.df$NREADS <- SAMPLE.SIZE[GENE.usage.df$SUBJECT]

  Binom.test.gene.cutoff <- Binom.test.gene.cutoff[[chain]]

  GENE.usage.df$min_frac <- sapply(1:nrow(GENE.usage.df), function(x){
    unique(Binom.test.gene.cutoff$min_frac[Binom.test.gene.cutoff$GENE==GENE.usage.df$GENE[x]])
  })

  GENE.usage.df.V <- GENE.usage.df %>% filter(grepl(paste0(chain,'V'),GENE))
  GENE.usage.df.J <- GENE.usage.df %>% filter(grepl(paste0(chain,'J'),GENE))

  GENE.usage.df.V$NREADS <- SAMPLE.SIZE.V[GENE.usage.df.V$SUBJECT]
  GENE.usage.df.J$NREADS <- SAMPLE.SIZE.J[GENE.usage.df.J$SUBJECT]

  GENE.usage.df.V$NREADS_SAMP <- SAMPLE.SIZE[GENE.usage.df.V$SUBJECT]
  GENE.usage.df.J$NREADS_SAMP <- SAMPLE.SIZE[GENE.usage.df.J$SUBJECT]

  GENE.loc.V <- GENE.loc[grep('V',GENE.loc)]
  GENE.usage.df.V <- binom_test_deletion(GENE.usage.df.V,cutoff=0.001,p.val.cutoff=0.01,chain=chain,GENE.loc.V)
  GENE.loc.J <- GENE.loc[grep('J',GENE.loc)]
  GENE.usage.df.J <- binom_test_deletion(GENE.usage.df.J,cutoff=0.005,p.val.cutoff=0.01,chain=chain,GENE.loc.J)

  if(chain=='IGH'){
    GENE.usage.df.D <- GENE.usage.df %>% filter(grepl(paste0(chain,'D'),GENE))
    GENE.usage.df.D$NREADS <- SAMPLE.SIZE.D[GENE.usage.df.D$SUBJECT]
    GENE.usage.df.D$NREADS_SAMP <- SAMPLE.SIZE[GENE.usage.df.D$SUBJECT]
    GENE.loc.D <- GENE.loc[grep('D',GENE.loc)]
    GENE.usage.df.D <- binom_test_deletion(GENE.usage.df.D,cutoff=0.005,p.val.cutoff=0.01,chain,GENE.loc.D)
    GENE.usage.df <- rbind(GENE.usage.df.V,GENE.usage.df.D,GENE.usage.df.J)
  }else{
    GENE.usage.df <- rbind(GENE.usage.df.V,GENE.usage.df.J)
  }
  GENE.usage.df <- GENE.usage.df %>% ungroup() %>% select(SUBJECT,GENE,FRAC,CUTOFF=min_frac,PVAL=pval_adj,DELETION=col)
  return(GENE.usage.df)
}

##########################################################################
#' Single chromosomal D or J gene deletions inferred by the V pooled method
#'
#' \code{deletionsByVpooled} detects D and J gene single chromosomal deletion.
#'
#' @details Function accepts a \code{data.frame} in Change-O format (\link{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{"SUBJECT"}: subject names
#'   \item \code{"V_CALL"}: V allele call(s) (in an IMGT format)
#'   \item \code{"D_CALL"}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{"J_CALL"}: J allele call(s) (in an IMGT format)
#' }
#'
#' @param  clip_db              a \code{data.frame} in Change-O format. See details.
#' @param  deletion_col         a vector of column names for which single chromosome deletions should be inferred. Default is J_CALL and D_CALL.
#' @param  count_thresh         integer, the minimun number of sequences mapped to a specific V gene to be included in the V pooled inference.
#' @param  deleted_genes         Double chromosome deletion summary table. A \code{data.frame} created by \code{binom_test_deletion}.
#' @param  min_minor_fraction   the minimum minor allele fraction to be used as an anchor gene. Default is 0.3
#'
#'
#' @return  data frame with single chromosome gene deletions
#'
#' @examples
#' # Load example data and germlines
#' data(sample_db)
#' # Infering V pooled deletions
#' del_db <- deletionsByVpooled(samples_db)
#' head(del_db)
#' @export
#not for light chain
deletionsByVpooled <- function(clip_db,deletion_col=c('D_CALL'),count_thresh=50,deleted_genes="",min_minor_fraction=0.3,
                               gene_sort='position',kThreshDel=3){


  if(!("SUBJECT" %in% names(clip_db))){clip_db$SUBJECT <- rep('S1',nrow(clip_db))}

  del.df <- c()
  for(sample_name in unique(clip_db$SUBJECT)){

    clip_db_sub <- clip_db[clip_db$SUBJECT==sample_name,]
    if(is.data.frame(deleted_genes)) deleted_genes_df <- deleted_genes %>% filter(SUBJECT==sample_name,grepl('IGHD|IGHJ',GENE))
    else deleted_genes_df <-c()
    ### Single V gene assignment
    #Number of iniial sequences
    nrows1 <- nrow(clip_db_sub)
    ### Only relevant genes with one assignment
    IND <- apply(clip_db_sub,1,function(x){sum(grep(',',x[c('V_CALL',deletion_col)],invert = F))})
    #IND <- names(IND[IND==0])
    clip_db_sub <- clip_db_sub[IND==0,]
    #Number of post multiple assignment reduction sequences
    nrows2 <- nrow(clip_db_sub)
    print(paste0(nrows1-nrows2,' sequnces were removed due to multiple assignments, ',nrows2,' sequences left'))

    ### Test for heterozygous V genes

    VGENES <- unique(sapply(strsplit(clip_db_sub$V_CALL,split = '*',fixed=T),'[',1))
    GENES <- unlist(sapply(VGENES,function(x){
      gene_counts <- table(grep(clip_db_sub$V_CALL,pattern = paste0(x,'*'),fixed = T,value = T));
      if(length(gene_counts) ==2 & (min(gene_counts)/sum(gene_counts)) >= min_minor_fraction & sum(gene_counts) >= count_thresh){return(x)}
    }))

    V.df <- list()
    if(length(GENES) > 0 ){
      print('The following genes used for pooled deletion detection')
      print(paste(GENES,sep=','))
      toHapGerm <- if(deletion_col=='D_CALL') HDGERM else HJGERM
      for(g in GENES){

        full.hap <- createFullHaplotype(clip_db_sub,toHap_col = deletion_col,hapBy_col = "V_CALL",hapBy = g,toHap_GERM = toHapGerm,
                                        relative_freq_priors = T,kThreshDel =kThreshDel,rmPseudo=T,
                                        deleted_genes = deleted_genes_df,chain='IGH',supress_print = T)

        full.hap$V_ALLELE_1 <- strsplit(names(full.hap)[5],'_')[[1]][2]
        full.hap$V_ALLELE_2 <- strsplit(names(full.hap)[6],'_')[[1]][2]

        names(full.hap)[5] <- paste0(g,'_',1)
        names(full.hap)[6] <- paste0(g,'_',2)
        V.df[[g]] <- rbind(V.df[[g]],full.hap)

      }
    } else {
      print('No heterozygous V genes found for deletion detection, try changing the parameters')
      exit()
    }

    df.compare <- c()

    GENES <- unlist(sapply( names(V.df),function(g){if(sample_name %in% V.df[[g]]$SUBJECT){return(g)}}))
    if(!is.null(GENES)){
      d.del.df <- c()

      for(g in GENES){
        tmp <- V.df[[g]] %>% filter(SUBJECT == sample_name,grepl('IGHD',GENE))
        tmp$DELETION <- apply(tmp,1,function(x){if(x[5]=='Unk' & x[6]!='Unk' | x[5]!='Unk' & x[6] == 'Unk'){return(1)};
          if(x[5]=='Del' & x[6]!='Del' | x[5]!='Del' & x[6] == 'Del'){return(2)};
          if(x[5]=='Unk' & x[6]=='Unk'){return(3)};
          if(x[5]=='Del' & x[6]=='Del'){return(4)};
          return(0)
        })
        tmp$K1_NEW <- sapply(1:nrow(tmp),function(i){if(!is.na(tmp$COUNTS2[i])){return(max(tmp$K1[i],tmp$K2[i]))} else {return(tmp$K1[i])}})
        tmp$COUNTS1_NEW <- sapply(1:nrow(tmp),function(i){if(!is.na(tmp$COUNTS2[i])){cnt1 <- as.numeric(strsplit(as.character(tmp$COUNTS1[i]),',')[[1]]);
        cnt2 <- as.numeric(strsplit(as.character(tmp$COUNTS2[i]),',')[[1]]);cnt <- cnt1+cnt2 ;  return( paste0(cnt,collapse=','))} else {return(as.character(tmp$COUNTS1[i]))}})

        tmp <- tmp %>% select(GENE,DELETION,K1_NEW,COUNTS1_NEW)
        tmp$V_GENE <- rep(g,nrow(tmp))
        tmp$DELETION2 <- ifelse(tmp$DELETION==0,0,1)
        d.del.df <- rbind(d.del.df,tmp)

      }


      tmp.df <- do.call('rbind',lapply(unique(d.del.df$GENE),function(x){
        x.df <- d.del.df %>% filter(GENE==x) %>%
          mutate(minCOUNT=(sapply(strsplit(as.character(COUNTS1_NEW),','),function(x) {min(as.numeric(x))})),
                 maxCOUNT=(sapply(strsplit(as.character(COUNTS1_NEW),','),function(x) {max(as.numeric(x))})));# %>%
        x.df <- x.df %>% filter(!is.na(K1_NEW));
        if(nrow(x.df) != 0){

          x.df <- x.df %>% group_by(DELETION2) %>% mutate(minCOUNTsum=sum(minCOUNT),maxCOUNTsum=sum(maxCOUNT)) %>% slice(1) %>%
            mutate(DELETION3=which.max(get_probabilites_with_priors(c(maxCOUNTsum,minCOUNTsum))[1:2]),
                   K2=max(get_probabilites_with_priors(c(maxCOUNTsum,minCOUNTsum))[1:2])-min(get_probabilites_with_priors(c(maxCOUNTsum,minCOUNTsum))[1:2])) ;
          return(x.df)
        }

      }))


      tmp.df.slct <- tmp.df %>% mutate(COUNTS2 = paste0(minCOUNTsum,',',maxCOUNTsum)) %>% select(GENE,DELETION2,DELETION3,K2,COUNTS2,V_GENE)
      tmp.df.slct$V_GENE <- ifelse(tmp.df.slct$DELETION2!=0,'POOLED_UNK','POOLED_KNOWN')
      tmp.df.slct <- tmp.df.slct[,-1]
      d.del.df <- d.del.df[names(d.del.df) != 'DELETION2']
      names(tmp.df.slct) <- names(d.del.df)
      tmp.df.slct$DELETION <- ifelse(tmp.df.slct$DELETION==2,0,1)
      ## ALL Vs
      tmp.df.slct.all <- tmp.df %>% ungroup() %>% group_by(GENE) %>% mutate(minCOUNTsum= sum(minCOUNTsum),maxCOUNTsum=sum(maxCOUNTsum)) %>% slice(1) %>% mutate(DELETION3=which.max(get_probabilites_with_priors(c(maxCOUNTsum,minCOUNTsum))[1:2]),K2=max(get_probabilites_with_priors(c(maxCOUNTsum,minCOUNTsum))[1:2])-min(get_probabilites_with_priors(c(maxCOUNTsum,minCOUNTsum))[1:2])) ;
      tmp.df.slct.all <- tmp.df.slct.all %>% mutate(COUNTS2 = paste0(minCOUNTsum,',',maxCOUNTsum)) %>% select(GENE,DELETION3,K2,COUNTS2,V_GENE)
      tmp.df.slct.all$V_GENE <- 'V(pooled)'
      names(tmp.df.slct.all) <- names(d.del.df)
      tmp.df.slct.all$DELETION <- ifelse(tmp.df.slct.all$DELETION==2,0,1)
      d.del.df <- rbind(d.del.df,tmp.df.slct)
      d.del.df <- rbind(d.del.df,as.data.frame(tmp.df.slct.all))
      d.del.df.pooled <-  d.del.df %>% filter(V_GENE=='V(pooled)')
      d.del.df.pooled$K1 <- round(as.numeric(d.del.df.pooled$K1),digits = 2)
      d.del.df.pooled$SUBJECT <- sample_name
      d.del.df.pooled$GENE <- factor(x = d.del.df.pooled$GENE,levels = GENE.loc[['IGH']])
     }
    del.df <- rbind(del.df,d.del.df.pooled)
  }
  return(del.df)

}
