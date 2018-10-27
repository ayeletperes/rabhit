# HaplotypeR graphic functions -----------------------------------------------------

#' @include rabhit.R
#' @include internal_functions.R
NULL

#' Graphical output of an inferred haplotype
#'
#' \code{plotHaplotype} visualizes an inferred haplotype.
#'
#' @details A \code{data.frame} in a haplotype format created by \code{createFullHaplotype} function.
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    html_output          If TRUE, a html5 interactive graph is outputed. Defualt is FALSE.
#' @param    gene_sort            If by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    text_size            the size of graph labels. Default is 14 (pts).
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#' @param    plotYaxis            if TRUE, Y axis labels (gene names) are plotted on the middle and right plots. Default is TRUE.
#' @param    chain                the Ig chain: IGH,IGK,IGL. Default is IGH.
#'
#'
#' @return   a haplotype visualization graph. If more than one subject is visualized, a pdf is created. If html_output is TRUE, a folder named html_output is created with individual graphs.
#'
#' @examples
#' # Load example data and germlines
#' data(sample_db, HVGERM, HDGERM)
#'
#' # Infering haplotype
#' hap_df = createFullHaplotype(sample_db,toHap_col=c("V_CALL","D_CALL"),hapBy_col="J_CALL",hapBy="IGHJ6",toHap_GERM=c(HVGERM,HDGERM));
#' plotHaplotype(hap_df)
#'
#' @export
plotHaplotype <- function(hap_table,html_output=FALSE, gene_sort = c("name", "position"),
                          text_size = 14, removeIGH=TRUE,plotYaxis=TRUE,chain=c('IGH','IGK','IGL'), ...)
{
  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(missing(gene_sort)) {
    gene_sort='position'
  }
  gene_sort <- match.arg(gene_sort)

  hapBy_cols = names(hap_table)[grep(chain,names(hap_table))]

  hapBy_alleles = gsub('_','*',hapBy_cols)

  if(!("SUBJECT" %in% names(hap_table))){hap_table$SUBJECT <- rep('S1',nrow(hap_table))}

  plot_list <- c()
  for(sample_name in unique(hap_table$SUBJECT)){

    GENE.loc.tmp <- GENE.loc[[chain]]
    hap_table_parse = parseHapTab(hap_table[hap_table$SUBJECT==sample_name,],hapBy_alleles)

    genotype = hap_table_parse$genotype
    kval.df = hap_table_parse$kval.df
    count.df = hap_table_parse$count.df

    alleles = strsplit(genotype$ALLELES, ",")

    geno2 = genotype
    r = 1
    for (g in 1:nrow(genotype)) {
      for (a in 1:length(alleles[[g]])) {
        geno2[r, ] = genotype[g, ]
        geno2[r, ]$ALLELES = alleles[[g]][a]
        r = r + 1
      }
    }
    if(gene_sort=='name'){
      geno2$GENE = factor(geno2$GENE, levels = rev(sortAlleles(unique(geno2$GENE), method = gene_sort)))
      kval.df$GENE = factor(kval.df$GENE, levels = rev(sortAlleles(unique(kval.df$GENE), method = gene_sort)))
    } else {

      names(GENE.loc.tmp) <- GENE.loc.tmp

      geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp))
      kval.df$GENE = factor(kval.df$GENE, levels = rev(GENE.loc.tmp))

      if(removeIGH){
        GENE.loc.tmp <- gsub('IG[H|K|L]','',GENE.loc.tmp)
        names(GENE.loc.tmp) <- GENE.loc.tmp
        geno2$GENE <- gsub('IG[H|K|L]','',geno2$GENE)
        kval.df$GENE <- gsub('IG[H|K|L]','',kval.df$GENE)
        geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp))
        kval.df$GENE = factor(kval.df$GENE, levels = rev(GENE.loc.tmp))
        geno2$hapBy <- gsub('IG[H|K|L]','',geno2$hapBy)
      } else {
        names(GENE.loc.tmp) <- GENE.loc.tmp

        geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp))
        kval.df$GENE = factor(kval.df$GENE, levels = rev(GENE.loc.tmp))

      }

    }

    geno2$hapBy <- gsub('_','*',geno2$hapBy,fixed = T)

    AlleleCol <- grep('[012]',unique(geno2$ALLELES),value = T,perl = T)
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

    AlleleCol<- names(c(alleles.comb,Unk='#dedede',Del='#6d6d6d'))
    names(AlleleCol) <- c(alleles.comb,Unk='#dedede',Del='#6d6d6d')

    transper <- sapply(AlleleCol,function(x){if(grepl('_',x)){mom_allele <- strsplit(x,'_')[[1]][1];
    all_novel <- grep(paste0(mom_allele,'_'),AlleleCol,value=T);
    if(length(all_novel)==1){return(0.5)};
    if(length(all_novel)==2){m=which(all_novel==x);return(ifelse(m==1,0.6,0.3))}
    if(length(all_novel)==3){m=which(all_novel==x);if(m==1){return(0.6)} ; return(ifelse(m==2,0.4,0.2))}
    } else (1)})
    names(transper) <- AlleleCol

    #remove 'mother' allele if added (when there is no germline allele but there is a novel)
    AlleleCol <- AlleleCol[AlleleCol %in% c(sort(grep('[012]',unique(geno2$ALLELES),value = T,perl = T)),'Unk','Del')]
    transper <- transper[names(transper) %in% AlleleCol ]

    p = ggplot(geno2, aes(x = GENE, fill = factor(ALLELES,levels=AlleleCol))) + theme_bw() +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            text = element_text(size = text_size), strip.background = element_blank(),
            strip.text = element_text(face = "bold"),axis.text = element_text(colour = "black"),
            panel.spacing = unit(0, "cm"),
            strip.switch.pad.grid = unit(0, "cm"),
            plot.margin = unit(c(0.25, 0, 0.2, 0), "cm")) + geom_bar(position = "fill",width = 0.9) +
      coord_flip() + xlab("") + ylab("") + facet_grid(paste0(".~", "hapBy"),switch = 'x') +
      scale_fill_manual(values=alpha(names(AlleleCol),transper),name='Alleles')

    if(!plotYaxis){
      p=p+theme(axis.text.y=element_blank())
    }



    kval.df$K_GROUPED <- bin_data(kval.df$K, bins=c(0, 1,2,3,4,5,10,20,50,Inf), binType = "explicit")
    ## plot K values
    pk <- ggplot(kval.df, aes(x = GENE, fill = K_GROUPED)) + theme_bw() +
      theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            text = element_text(size = text_size), strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(0, "cm"),
            strip.switch.pad.grid = unit(0, "cm"),
            plot.margin = unit(c(0.25, 0, 0.2, 0), "cm")) + geom_bar(position = "fill",width = 0.7)  +
      coord_flip() +xlab("")+ ylab("") + facet_grid(paste0(".~", "hapBy"),switch = 'x')


    if(!html_output){
      p.legend <- get_legend(p)
      p = p +theme(legend.position='none')

      pk = pk + scale_fill_brewer(name=expression('log'[10]*'(lK)'),drop = FALSE)
      pk.legend <- get_legend(pk)
      pk = pk +theme(legend.position='none')

      p = p + do.call(theme, list(...)) + theme( axis.title.x=element_blank())
      pk = pk + do.call(theme, list(...)) + theme( axis.title=element_blank())

      p.legends  <-  plot_grid(pk.legend,p.legend, ncol=1, rel_heights = c(0.2,1),align = 'hv')


      p1 <- ggdraw(plot_grid(p,pk,p.legends, ncol=3, rel_widths=c(0.1,0.1 ,0.1)))
    }else{
      p = p + do.call(theme, list(...)) + theme( axis.title.x=element_blank())
      pk = pk + do.call(theme, list(...)) + theme( axis.title=element_blank())
    }



    ########################################################################################################

    ## Plot coutns
    if(gene_sort=='name'){
      count.df$GENE = factor(count.df$GENE, levels = rev(sortAlleles(unique(count.df$GENE), method = gene_sort)))

      if(removeIGH){
        count.df$GENE <- gsub('IG[H|K|L]','',count.df$GENE)
        count.df$GENE = factor(count.df$GENE, levels = rev(sortAlleles(unique(count.df$GENE), method = gene_sort)))

      }
    } else {
      if(removeIGH){
        GENE.loc.tmp <- gsub('IG[H|K|L]','',GENE.loc.tmp)
        names(GENE.loc.tmp) <- GENE.loc.tmp

        count.df$GENE <- gsub('IG[H|K|L]','',count.df$GENE)
        count.df$GENE = factor(count.df$GENE, levels = rev(GENE.loc.tmp))

      } else {
        names(GENE.loc.tmp) <- GENE.loc.tmp

        count.df$GENE = factor(count.df$GENE, levels = rev(GENE.loc.tmp))
      }
    }

    ## TO visualy make coutns of 1 not look like 0 , one is added

    count.df$COUNT2 <- ifelse(count.df$HapBy == hapBy_alleles[1],-1*log10(as.numeric(count.df$COUNT)+1),log10(as.numeric(count.df$COUNT)+1))
    count.df$COUNT2[count.df$COUNT2==Inf |count.df$COUNT2== -Inf  ] <- 0


    AlleleCol <- c(sort(grep('[012]',unique(count.df$ALLELES),value = T,perl = T)),'Unk','Del')


    p2 <- ggplot(count.df,aes(x=GENE,y=COUNT2,fill=factor(ALLELES,levels=AlleleCol))) + geom_bar(stat="identity",position = 'Dodge',width = 0.9) +
      coord_flip() + background_grid(minor='none') +
      theme(legend.position="none",strip.text = element_text(face = "bold"),axis.text = element_text(colour = "black"),
            text = element_text(size = text_size),plot.margin = unit(c(0.25, 0, -0.05, 0), "cm")) + scale_y_continuous(breaks = seq(-3, 3, by = 1),labels = c(3:0,1:3)) +
      ylab(expression('log'[10]*'(Count+1)')) + xlab('Gene')  + geom_hline(yintercept=c(0), linetype="dotted")

    if(is.null(ALLELE_PALETTE)){
      AlleleCol <- c(sort(grep('[012]',unique(count.df$ALLELES),value = T,perl = T)),'Unk','Del')
      tmp.col <- scale_fill_hue(h.start = 10, h=c(10, 270))
      len.col <- length(AlleleCol)-2
      names(AlleleCol) <-  c(tmp.col$palette(len.col),'#dedede','#6d6d6d')
      p2 = p2 +scale_fill_manual(values=(names(AlleleCol)),name='ALLELES')
    } else {
      AlleleCol <- grep('[012]',unique(count.df$ALLELES),value = T,perl = T)
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


      AlleleCol<- names(c(alleles.comb,Unk='#dedede',Del='#6d6d6d'))
      names(AlleleCol) <- c(alleles.comb,Unk='#dedede',Del='#6d6d6d')

      transper <- sapply(AlleleCol,function(x){if(grepl('_',x)){mom_allele <- strsplit(x,'_')[[1]][1];
      all_novel <- grep(paste0(mom_allele,'_'),AlleleCol,value=T);
      if(length(all_novel)==1){return(0.5)};
      if(length(all_novel)==2){m=which(all_novel==x);return(ifelse(m==1,0.6,0.3))}
      if(length(all_novel)==3){m=which(all_novel==x);if(m==1){return(0.6)} ; return(ifelse(m==2,0.4,0.2))}
      } else (1)})
      names(transper) <- AlleleCol

      #remove 'mother' allele if added (when there is no germline allele but there is a novel)
      AlleleCol <- AlleleCol[AlleleCol %in% c(sort(grep('[012]',unique(count.df$ALLELES),value = T,perl = T)),'Unk','Del')]
      transper <- transper[names(transper) %in% AlleleCol ]

      p2 = p2 +scale_fill_manual(values=alpha(names(AlleleCol),transper),name='ALLELES')

    }

    ### Plot both panels

    if(html_output){
      pk <- pk +   scale_fill_brewer(name='log<sub>10</sub>(lK)',drop = FALSE)
      p2 <- p2 + ylab('log<sub>10</sub>(Count+1)')
      p.l <- ggplotly(p,height = 1000,width = 700) %>% plotly::layout(showlegend=FALSE)
      pk.l <- ggplotly(pk,height = 1000,width = 700) %>% plotly::layout(showlegend=TRUE)

      pk.l$x$layout$annotations[[1]]$text = p.l$x$layout$annotations[[1]]$text
      pk.l$x$layout$annotations[[2]]$text = p.l$x$layout$annotations[[2]]$text

      p2.l <- ggplotly(p2,height = 1300,width = 1000) %>% plotly::layout(margin=list(b=50),
                                                                 yaxis = list(title = paste0(c(rep("&nbsp;", 3),
                                                                                               "Gene",
                                                                                               rep("&nbsp;", 3),
                                                                                               rep("\n&nbsp;", 1)),
                                                                                             collapse = "")),
                                                                 showlegend=TRUE)

      p2.l$x$layout$xaxis$ticktext = c(lapply(p2.l$x$layout$xaxis$ticktext[1:match('0',p2.l$x$layout$xaxis$ticktext)-1],function(x) paste0('-',x)),
                                       p2.l$x$layout$xaxis$ticktext[match('0',p2.l$x$layout$xaxis$ticktext):length(p2.l$x$layout$xaxis$ticktext)])

      p2.l$x$data[[length(p2.l$x$data)]]$y[2] = p2.l$x$layout$yaxis$range[2]

      mgsub <- function(pattern, replacement, x, ...) {
        if (length(pattern)!=length(replacement)) {
          stop("pattern and replacement do not have the same length.")
        }
        result <- x
        for (i in 1:length(pattern)) {
          result <- gsub(pattern[i], replacement[i], result, ...)
        }
        result
      }


      text_for_hovertext <- function(labels,count.df) {
        for(i in 1:length(labels)){
          label <- labels[i]
          gene <- strsplit(strsplit(label,'<')[[1]][1],' ')[[1]][2]
          allele <- strsplit(label,'Allele: ')[[1]][2]
          if(!is.na(allele)){
            count <- as.numeric(strsplit(strsplit(label,'<br />Count: ')[[1]][2],'<')[[1]][1])
            if(count%%1!=0) count <- count.df %>% filter(GENE==gene&ALLELES==allele&round(COUNT3,nchar(as.character(count))-2)==count) %>% select(COUNT)
            else count <- count.df %>% filter(GENE==gene&ALLELES==allele&COUNT3==count) %>% select(COUNT)
            labels[i] <- paste0('Gene: ',gene,'<br />Allele: ',allele,'<br />Count: ',count[1,])}
        }
        return(labels)
      }

      count.df$COUNT3 <- abs(count.df$COUNT2)

      for(i in 1:length(p2.l$x$data)){
        p2.l$x$data[[i]]$text <- mgsub(c("~GENE","~COUNT2","~factor[(]ALLELES[,] levels [=] AlleleCol[)]","~yintercept: 0"),
                                       c("Gene","Count","Allele",""),p2.l$x$data[[i]]$text)
        p2.l$x$data[[i]]$text <- text_for_hovertext(p2.l$x$data[[i]]$text,count.df)
      }

      for(i in 1:length(p.l$x$data)){
        p.l$x$data[[i]]$text <- mgsub(c("~GENE","~factor[(]ALLELES[,] levels [=] AlleleCol[)]"),
                                      c("Gene","Allele"),p.l$x$data[[i]]$text)
      }

      for(i in 1:length(pk.l$x$data)){
        pk.l$x$data[[i]]$text <- mgsub(c("~GENE","~K[_]GROUPED"),
                                       c("Gene","K"),pk.l$x$data[[i]]$text)
      }


      p.l.c <- suppressWarnings(subplot(p2.l,p.l,pk.l,widths = c(0.4,0.2,0.2),shareY = T,titleX = TRUE,margin = 0.01,which_layout = 1))
      p.l.c$x$layout$annotations[[6]]$text = "log<sub>10</sub>(lK)-------------"
      p.l.c$x$layout$annotations[[6]]$xanchor="center"
      p.l.c$x$layout$annotations[[6]]$y=0.99-0.0234*(length(AlleleCol)+2.5)#0.52
      p.l.c$x$layout$annotations[[6]]$x=1.02


      p.l.c$x$layout$annotations[[3]] <- p.l.c$x$layout$annotations[[6]]
      p.l.c$x$layout$annotations[[3]]$text = "Alleles-------------"
      p.l.c$x$layout$annotations[[3]]$y=0.99
      p.l.c$x$layout$annotations[[3]]$x=1.02
      p.l.c$x$layout$annotations[[3]]$legendTitle=FALSE

      for(i in (length(AlleleCol)-2+1):(length(p.l.c$x$data)-(length(pk.l$x$data)+4))){
        p.l.c$x$data[[i]]$showlegend <- FALSE
      }

      plot_list[[sample_name]] <- p.l.c

    }else{

    p <- plot_grid(p2+ theme(plot.margin = unit(c(0.5,0,0,0),"lines")), p1+theme(plot.margin = unit(c(0,0,0.5,0),"lines")),nrow=1,rel_widths = c(1,2))
    # now add the title
    title <- ggdraw() + draw_label(sample_name, fontface='bold')

    plot_list[[sample_name]] <- plot_grid(title, p, ncol=1, rel_heights=c(0.05, 1))}

  }
  if(length(plot_list) != 1){

    if(html_output){
      dir.create(file.path(getwd(),'html_output'))
      for(sample_name in names(plot_list)){
        htmlwidgets::saveWidget(plot_list[[sample_name]], paste0(getwd(),'/html_output/',sample_name,".html"),selfcontained = T)
      }
    }
    else{
      pdf(paste0(getwd(),'/haplotype_output.pdf'),height = 20,width = 15)
      for(p in plot_list){
        plot(p)
      }
      dev.off()
    }

  }
  else if(html_output) return(plot_list[[1]]) else plot(plot_list[[1]])
}
########################################################################################################
#' Graphical output of single chromosome deletions
#'
#' \code{deletionHeatmap} creates a graphical output of the single chromosome deletions in multiple samples.
#'
#' @details A \code{data.frame} created by \code{createFullHaplotype}.
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    html_output          If TRUE, a html5 interactive graph is outputed insteaed of the normal plot. Defualt is FALSE
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#'
#'
#' @return   a single chromosome deletion visualization.
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db, HVGERM, HDGERM)
#'
#' # Infering haplotype
#' hap_df = createFullHaplotype(samples_db,toHap_col=c("V_CALL","D_CALL"),hapBy_col="J_CALL",hapBy="IGHJ6",toHap_GERM=c(HVGERM,HDGERM));
#' deletionHeatmap(hap_df)
#'
#' @export
deletionHeatmap <- function(hap_table,html_output=FALSE,chain=c('IGH','IGK','IGL'), ...){

  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(!("SUBJECT" %in% names(hap_table))){hap_table$SUBJECT <- rep('S1',nrow(hap_table))}

  hapBy_cols = names(hap_table)[grep(chain,names(hap_table))]

  genes_hap <- unique(substr(hap_table$GENE,4,4))

  GENE.loc <- GENE.loc[[chain]][grep(paste0(genes_hap,collapse = '|'),GENE.loc[[chain]])]

  ALLELE_01_col = hapBy_cols[1]
  ALLELE_02_col = hapBy_cols[2]

  ALLELE_01_num = strsplit(ALLELE_01_col,'_')[[1]][2]
  ALLELE_02_num = strsplit(ALLELE_02_col,'_')[[1]][2]

  ### create deletion streches heatmap
  hap_table$K1[is.na(hap_table$K1)] <- 0
  hap_table$K2[is.na(hap_table$K2)] <- 0

  hap_table.del.heatmap <- hap_table %>% rowwise %>%
    mutate(K=max(as.numeric(K1),as.numeric(K2),na.rm = T)) %>% select_(.dots=c('SUBJECT','GENE',ALLELE_01_col,'K'))

  hap_table.del.heatmap$HapBy <- rep(ALLELE_01_num,nrow(hap_table.del.heatmap))
  names(hap_table.del.heatmap)[3] <- ALLELE_02_col

  hap_table.del.heatmap <- rbind(hap_table.del.heatmap,data.frame(hap_table %>% rowwise %>%
                                                                    mutate(K=max(as.numeric(K1),as.numeric(K2),na.rm = T)) %>%
                                                                    select_(.dots=c('SUBJECT','GENE',ALLELE_02_col,'K')),HapBy=ALLELE_02_num))

  names(hap_table.del.heatmap)[3] <- 'ALLELE'

  hap_table.del.heatmap$K[hap_table.del.heatmap$K == Inf] <- 0
  hap_table.del.heatmap$K[hap_table.del.heatmap$K == -Inf] <- 0
  hap_table.del.heatmap$DEL <- ifelse(hap_table.del.heatmap$ALLELE == 'Del'   ,2,0 )
  hap_table.del.heatmap$DEL[hap_table.del.heatmap$ALLELE != 'Del' & hap_table.del.heatmap$K<=3] <- 4
  hap_table.del.heatmap$DEL[hap_table.del.heatmap$ALLELE == 'Unk'] <- 1
  hap_table.del.heatmap$DEL[hap_table.del.heatmap$ALLELE == 'NA'] <- 3


  #manual reshape
  hap_table.del.heatmap.02 <- hap_table.del.heatmap[hap_table.del.heatmap$HapBy==ALLELE_01_num,]
  hap_table.del.heatmap.03 <- hap_table.del.heatmap[hap_table.del.heatmap$HapBy==ALLELE_02_num,]
  hap_table.del.heatmap.02$GENE2 <- gsub('IG[H|K|L]','',hap_table.del.heatmap.02$GENE)
  hap_table.del.heatmap.03$GENE2 <- gsub('IG[H|K|L]','',hap_table.del.heatmap.03$GENE)

  hap_table.del.heatmap.02$SUBJECT <- factor(x = hap_table.del.heatmap.02$SUBJECT,
                                             levels = unique(hap_table.del.heatmap.02$SUBJECT))

  hap_table.del.heatmap.02$GENE2 <- factor(x = hap_table.del.heatmap.02$GENE2,
                                           levels = gsub('IG[H|K|L]','',GENE.loc))

  hap_table.del.heatmap.03$SUBJECT <- factor(x = hap_table.del.heatmap.03$SUBJECT,
                                             levels = unique(hap_table.del.heatmap.03$SUBJECT))

  hap_table.del.heatmap.03$GENE2 <- factor(x = hap_table.del.heatmap.03$GENE2,
                                           levels = gsub('IG[H|K|L]','',GENE.loc))
  heatmap.df <- rbind(hap_table.del.heatmap.02,hap_table.del.heatmap.03)

  ALLELE_01_col = gsub('_','*',gsub('IG[H|K|L]','',ALLELE_01_col))
  ALLELE_02_col = gsub('_','*',gsub('IG[H|K|L]','',ALLELE_02_col))

  heatmap.df$HapBy <- ifelse(heatmap.df$HapBy==ALLELE_01_num,ALLELE_01_col,ALLELE_02_col)
  heatmap.plot <- ggplot(data = heatmap.df, aes(x = GENE2, y = SUBJECT)) + theme_bw()+
    geom_tile(aes(fill = as.character(DEL))) + facet_wrap(~HapBy,nrow=2)+ scale_x_discrete(drop=FALSE)+
    scale_fill_manual(name='lK',labels=c('0','1','2','3','4'),values = c('white','lightblue','blue','grey40','grey90'))  + ylab('Subject') + xlab('Gene') +
    theme(strip.text = element_text(size=18),axis.title = element_text(size=18),axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90,vjust = 0.5 ,hjust = 1),plot.margin = margin(b = 30),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom")

  del.df.heatmap <- heatmap.df
  del.df.heatmap <- del.df.heatmap %>% filter(ALLELE=='Del')

  del.df.heatmap.cnt <- del.df.heatmap %>% group_by(SUBJECT,GENE2) %>% mutate(n=n()) %>% dplyr::slice(1)
  del.df.heatmap.cnt$HapBy[del.df.heatmap.cnt$n==2] <- "Both"
  del.df.heatmap.cnt <- del.df.heatmap.cnt %>% group_by(GENE2,HapBy) %>% count_()

  del.df.heatmap.cnt$HapBy <- factor(del.df.heatmap.cnt$HapBy,levels=c(ALLELE_01_col,ALLELE_02_col,'Both'))
  pdel <- ggplot(del.df.heatmap.cnt,aes(x=GENE2,y=nn,fill=(HapBy))) + theme_bw() +
    geom_bar(stat='identity',position = 'stack')  +
    theme(strip.background = element_blank(),
          axis.text = element_text(size=14),
          axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),plot.margin = margin(0,8,0,7,"pt"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    ylab('Number of individuals\nwith a deletion')+
    scale_fill_manual(name = "Chromosome",values = c('darksalmon','deepskyblue','darkolivegreen3','grey50'))  + scale_x_discrete(drop=FALSE) + xlab('')

  if(html_output){
    m = list(l = 200,r = 40,b = 150,t = 50,pad = 0)

    pdel.l <- ggplotly(pdel,height = 900,width = 1500) %>% plotly::layout(autosize = F, margin = m)


    heatmap.plot.l <- ggplotly(heatmap.plot,height = 900,width = 1500) %>% plotly::layout(autosize = F, margin = m,xaxis=list(showgrid=T))

    heatmap.bar.l <- subplot(pdel.l,heatmap.plot.l,nrows = 2) %>%
      plotly::layout(xaxis = list(domain=list(x=c(0,0.5),y=c(0,0.5))),
             xaxis2 = list(domain=list(x=c(0.5,1),y=c(0.5,1))))

  }
  legend <- cowplot::get_legend(heatmap.plot)
  heatmap.plot <- heatmap.plot + theme(legend.position = 'none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
  pdel <- pdel + theme(legend.position = c(0.9, 0.8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))



  comb <- plot_grid(pdel, heatmap.plot,ncol=1,rel_heights=c(0.15, 0.3),align = "h")
  plot(plot_grid(comb, legend,nrow=2,rel_heights=c(1, 0.1)))
}

########################################################################################################
#' Graphical output of alleles division by chromosome
#'
#' \code{hapHeatmap} creates a graphical output of the alleles per gene in chromosome in multiple samples.
#'
#' @details A \code{data.frame} created by \code{createFullHaplotype}.
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    gene_sort            If by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#'
#' @return   a single chromosome deletion visualization.
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db, HVGERM, HDGERM)
#'
#' # Infering haplotype
#' hap_df = createFullHaplotype(samples_db,toHap_col=c("V_CALL","D_CALL"),hapBy_col="J_CALL",hapBy="IGHJ6",toHap_GERM=c(HVGERM,HDGERM));
#' hapHeatmap(hap_df)
#'
#' @export
hapHeatmap <- function(hap_table,chain=c('IGH','IGK','IGL'),gene_sort='position',removeIGH=TRUE, ...){

  hapBy_alleles <-  gsub('_','*',names(hap_table)[grep(chain,names(hap_table))])
  hapBy_cols <- gsub('IG[H|K|L]','',hapBy_alleles)

  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  genos <- c()
  for(sample_name in unique(hap_table$SUBJECT)){

    hap_table_parse = parseHapTab(hap_table[hap_table$SUBJECT==sample_name,],hapBy_alleles)
    GENE.loc.tmp <- GENE.loc[[chain]]

    genotype = hap_table_parse$genotype
    kval.df = hap_table_parse$kval.df
    count.df = hap_table_parse$count.df

    alleles = strsplit(genotype$ALLELES, ",")

    geno2 = genotype
    r = 1
    for (g in 1:nrow(genotype)) {
      for (a in 1:length(alleles[[g]])) {
        geno2[r, ] = genotype[g, ]
        geno2[r, ]$ALLELES = alleles[[g]][a]
        r = r + 1
      }
    }
    if(gene_sort=='name'){
      geno2$GENE = factor(geno2$GENE, levels = rev(sortAlleles(unique(geno2$GENE), method = gene_sort)))
    } else {

      names(GENE.loc.tmp) <- GENE.loc.tmp

      geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp))

      if(removeIGH){
        GENE.loc.tmp <- gsub('IG[H|K|L]','',GENE.loc.tmp)
        names(GENE.loc.tmp) <- GENE.loc.tmp
        geno2$GENE <- gsub('IG[H|K|L]','',geno2$GENE)
        geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp))
        geno2$hapBy <- gsub('IG[H|K|L]','',geno2$hapBy)
      } else {
        names(GENE.loc.tmp) <- GENE.loc.tmp

        geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp))

      }

    }

    geno2$hapBy <- gsub('_','*',geno2$hapBy,fixed = T)
    geno2$SUBJECT <- sample_name
    genos <- rbind(genos,geno2)
  }

  heatmap.df <- genos %>% group_by(SUBJECT,hapBy,GENE) %>% mutate(n=n())
  heatmap.df$freq <-ifelse(heatmap.df$n==2,0.5,1)
  heatmap.df$GENE <- factor(heatmap.df$GENE, levels =  gsub('IG[H|K|L]','',GENE.loc[[chain]]))
  AlleleCol <- grep('[012]',unique(heatmap.df$ALLELES),value = T,perl = T)
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

  AlleleCol<- names(c(alleles.comb,Unk='#dedede',Del='#6d6d6d'))
  names(AlleleCol) <- c(alleles.comb,Unk='#dedede',Del='#6d6d6d')

  transper <- sapply(AlleleCol,function(x){
    if(grepl('_',x)){
      mom_allele <- strsplit(x,'_')[[1]][1];
      all_novel <- grep(paste0(mom_allele,'_'),AlleleCol,value=T);
      if(length(all_novel)==1){return(0.5)};
      if(length(all_novel)==2){m=which(all_novel==x);return(ifelse(m==1,0.6,0.3))}
      if(length(all_novel)==3){m=which(all_novel==x);if(m==1){return(0.6)} ; return(ifelse(m==2,0.4,0.2))}
      if(length(all_novel)>3){m=which(all_novel==x);if(m==1){return(0.85)} ; return(0.85-m/10)}
    }else (1)})
  names(transper) <- AlleleCol

  #remove 'mother' allele if added (when there is no germline allele but there is a novel)
  AlleleCol <- AlleleCol[AlleleCol %in% c(sort(grep('[012]',unique(heatmap.df$ALLELES),value = T,perl = T)),'Unk','Del')]
  AlleleCol.tmp <- names(AlleleCol)
  names(AlleleCol.tmp) <- AlleleCol
  transper <- transper[names(transper) %in% AlleleCol ]

  heatmap.df$title <- ifelse(heatmap.df$hapBy==hapBy_cols[1],hapBy_cols[1],hapBy_cols[2])

  col <- ifelse(c(1:(nrow(unique(heatmap.df[heatmap.df$hapBy==hapBy_cols[1],'SUBJECT']))*4))%%3==0,'black','white')

  p <- ggplot(heatmap.df[heatmap.df$hapBy==hapBy_cols[1],], (aes(x = GENE, y = freq, fill = factor(ALLELES,levels=AlleleCol)))) +
    geom_col(position = "fill", width = 0.95) +
    scale_fill_manual(values=alpha(names(AlleleCol),transper),name='Alleles',drop=F)+
    facet_grid(SUBJECT ~ title, as.table = FALSE, switch = "y") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +  theme(axis.text.x = element_text(angle = 90,vjust = 0.5 ,hjust = 1, size=10,colour = 'black')
                                                , strip.text.x = element_text(size=14)
                                                , strip.text.y = element_text(angle = 180, size=10)
                                                , panel.grid = element_blank()
                                                , strip.placement = "outside"
                                                , axis.ticks.y = element_line(colour = col)
                                                , axis.text.y = element_blank()
                                                , strip.background.y = element_blank()
                                                , panel.spacing.y = unit(0.9, "pt")) + labs(y='',x='')

  col <- ifelse(c(1:(nrow(unique(heatmap.df[heatmap.df$hapBy==hapBy_cols[2],'SUBJECT']))*4))%%3==0,'black','white')

  p1 <- ggplot(heatmap.df[heatmap.df$hapBy==hapBy_cols[2],], (aes(x = GENE, y = n, fill = factor(ALLELES,levels=AlleleCol)))) +
    geom_col(position = "fill", width = 0.95) +
    scale_fill_manual(values=alpha(names(AlleleCol),transper),name='Alleles',drop=F)+
    facet_grid(SUBJECT ~ title, as.table = FALSE, switch = "y") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +  theme(axis.text.x = element_text(angle = 90,vjust = 0.5 ,hjust = 1, size=10,colour = 'black')
                                                , strip.text.x = element_text(size=14)
                                                , strip.text.y = element_text(angle = 180, size=10)
                                                , panel.grid = element_blank()
                                                , strip.placement = "outside"
                                                , axis.ticks.y = element_line(colour = col)
                                                , axis.text.y = element_blank()
                                                , strip.background.y = element_blank()
                                                , panel.spacing.y = unit(0.9, "pt")) + labs(y='')

  legend <- get_legend(p)
  plot(plot_grid(plot_grid(p+theme(legend.position = 'none'),
                                        p1+theme(legend.position = 'none'),
                                        nrow=2,align = 'hv'),legend,rel_widths = c(0.7,0.15),ncol=2))
}
########################################################################################################
#' Graphical output of double chromosome deletions
#'
#' \code{plotDeletionsByBinom} creates a graphical output of the double chromosome deletions in multiple samples.
#'
#' @details A \code{data.frame} created by \code{binom_test_deletion}.
#'
#' @param    GENE.usage.df        Double chromosome deletion summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    genes.low.cer        a vector of IGH genes known to be with low certantiny in the binomial test. Default is IGHV3-43 and IGHV3-20
#' @param    genes.dup            a vector of IGH genes known to have a duplicated gene. Default is IGHD4-11 that his duplicate is IGHD4-4 and IGHD5-18 that his duplicate is IGHD5-5
#'
#' @return   a double chromosome deletion visualization.
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db)
#'
#' # Infering haplotype
#' hap_df = binom_test_deletion(samples_db);
#' plotDeletionsByBinom(hap_df)
#'
#' @export
plotDeletionsByBinom <- function(GENE.usage.df,chain=c('IGH','IGK','IGL'),genes.low.cer=c('IGHV3-43','IGHV3-20'), genes.dup=c('IGHD4-11','IGHD5-18'),...){

  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(!("SUBJECT" %in% names(GENE.usage.df))){GENE.usage.df$SUBJECT <- rep('S1',nrow(GENE.usage.df))}

  genes_hap <- unique(substr(GENE.usage.df$GENE,4,4))
  GENE.loc <- GENE.loc[[chain]][GENE.loc[[chain]] %in% GENE.usage.df$GENE]

  GENE.usage.df$GENE2 <- factor(gsub(chain,'',GENE.usage.df$GENE), levels=gsub(chain,'',GENE.loc))

  colvec <- ifelse(GENE.loc%in%genes.low.cer, "red", ifelse(GENE.loc%in%genes.dup, "purple", "black"))

  ### gene usage with deletions in population according to binom test
  p.del <- ggplot(GENE.usage.df,aes(x=GENE2,y=FRAC)) + geom_boxplot(outlier.colour=NA) +
    geom_jitter(aes(x=GENE2,color=DELETION),width = 0.25,size=0.5)+ theme(axis.text.y = element_text(size=16),
                                                                     axis.title = element_text(size=16),
                                                                     axis.text.x = element_text(size=14,angle = 90, hjust = 1,vjust=0.5,color=colvec),
                                                                     legend.text = element_text(size=16),
                                                                     legend.position = 'none')+
    ylab('Fraction') + xlab('')  + scale_color_manual(name='',labels=c('Deletion','No Deletion','NA'),values = c('blue','black','grey40'),drop=F)+
    guides(color = guide_legend(override.aes = list(size=5)))

  ### heat map of deletions in population according to binom test
  GENE.usage.df$DELETION <- factor(GENE.usage.df$DELETION,levels=levels(GENE.usage.df$DELETION))

  heatmap.plot <-   ggplot(data = GENE.usage.df, aes(x = GENE2, y = SUBJECT)) +
    geom_tile(aes(fill = DELETION)) +
    scale_fill_manual(name='',labels=c('Deletion','No Deletion','NA'),values = c('blue','white','grey40'),drop=F)+
    scale_x_discrete(drop=FALSE) +ylab('Subject')+ xlab('')+
    theme(axis.text = element_text(size = 12),axis.title.y = element_text(size=18),
          axis.text.x = element_text(angle = 90,vjust = 0.5 ,hjust = 1,size=14),
          legend.key = element_rect(colour = 'black' , size = 0.5, linetype='solid'),legend.text = element_text(size=16),
          legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom")

  legend <- cowplot::get_legend(heatmap.plot)
  heatmap.plot <- heatmap.plot + theme(legend.position = 'none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p.del <- p.del + theme(legend.position = 'none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"))

  comb <- plot_grid(p.del, heatmap.plot,ncol=1,rel_heights=c(0.15, 0.3),align = "hv")
  plot(plot_grid(comb, legend,nrow=2,rel_heights=c(1, 0.2)))
}

##########################################################################
#' Graphical output for single chromosome D or J gene deletions according to V pooled method
#'
#' \code{plotDeletionsByVpooled} graphical output for single chromosome D or J gene deletions (for heavy chain only).
#'
#' @details A \code{data.frame} created by \code{deletionsByVpooled}.
#'
#' @param  del.df  A \code{data.frame} created by \code{deletionsByVpooled}.
#' @param  K_ranges vector of one or two integers for log(K) certainty level thresholds
#'
#' @return   a single chromosome deletion visualization.
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db)
#' del_db <- deletionsByVpooled(samples_db)
#' plotDeletionsByVpooled(del_db)
#' @export


plotDeletionsByVpooled <- function(del.df,K_ranges=c(3,7)){


  if(!("SUBJECT" %in% names(del.df))){del.df$SUBJECT <- rep('S1',nrow(del.df))}

  if(length(K_ranges)==2){
    del.df$EVENT <- unlist(sapply(1:nrow(del.df), function(i){
      if(del.df$DELETION[i] > 0 & del.df$K1[i] < K_ranges[1]) return(1)
      if(del.df$DELETION[i] > 0  &  del.df$K1[i] >= K_ranges[1] &del.df$K1[i] < K_ranges[2]) return(2)
      if(del.df$DELETION[i] > 0  &  del.df$K1[i] > K_ranges[2]) return(3)
      if(del.df$DELETION[i] == 0 & del.df$K1[i] < K_ranges[1]) return(4)
      if(del.df$DELETION[i] == 0  &  del.df$K1[i] >= K_ranges[1] &del.df$K1[i] < K_ranges[2]) return(5)
      if(del.df$DELETION[i] == 0  &  del.df$K1[i] > K_ranges[2]) return(6)
    }))


    del.df$EVENT <- factor(del.df$EVENT,levels=1:6)
    del.df$GENE2 <- factor(gsub('IGH','',del.df$GENE),levels = gsub('IGH','',GENE.loc[['IGH']]))
    labels1=c(paste0('Deletion lK<',K_ranges[1]),paste0('Deletion ',K_ranges[1],'<=lK<',K_ranges[2]),
              paste0('Deletion lK>=',K_ranges[2]),paste0('No deletion lK<',K_ranges[1]),paste0('No deletion ',K_ranges[1],'<=lK<',K_ranges[2]),
              paste0('No deletion lK>=',K_ranges[2]))
    values1 = c('lightblue','cornflowerblue','black','lightpink','lightcoral','lightgrey')
  }

  if(length(K_ranges)==1){
    del.df$EVENT <- unlist(sapply(1:nrow(del.df), function(i){
      if(del.df$DELETION[i] > 0 & del.df$K1[i] < K_ranges[1]) return(1)
      if(del.df$DELETION[i] > 0  &  del.df$K1[i] >= K_ranges[1]) return(2)
      if(del.df$DELETION[i] == 0 & del.df$K1[i] < K_ranges[1]) return(3)
      if(del.df$DELETION[i] == 0  &  del.df$K1[i] >= K_ranges[1]) return(4)
    }))


    del.df$EVENT <- factor(del.df$EVENT,levels=1:4)
    del.df$GENE2 <- factor(gsub('IGH','',del.df$GENE),levels = gsub('IGH','',GENE.loc[['IGH']]))
    labels1=c(paste0('Deletion lK<',K_ranges[1]),
              paste0('Deletion lK>=',K_ranges[1]),paste0('No deletion lK<',K_ranges[1]),
              paste0('No deletion lK>=',K_ranges[1]))
    values1 = c('lightblue','black','lightpink','lightgrey')
  }

  heatmap.plot <- ggplot(data = del.df, aes(x = GENE2, y = SUBJECT)) +
    geom_tile(aes(fill = EVENT)) +
    scale_fill_manual(name='',labels=labels1,
                      values = values1,drop=F) + ylab('Subject') + xlab('Gene') +
    theme(strip.text = element_text(size=14),axis.title = element_text(size=14),axis.text = element_text(size = 12),
          legend.position = 'bottom',axis.text.x = element_text(angle = 90,vjust = 0.5 ,hjust = 1),legend.text=element_text(size=14))



  plot(heatmap.plot)

}
