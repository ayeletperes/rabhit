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
plotHaplotype <- function(hap_table, html_output=FALSE, gene_sort = c("name", "position"),
                          text_size = 14, removeIGH=TRUE, plotYaxis=TRUE, chain = c('IGH','IGK','IGL'), ...)
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

    haplo.db <- parseHapTab(hap_table[hap_table$SUBJECT==sample_name,], chain = chain)
    geno.df <- sortDFByGene(haplo.db$geno.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    kval.df <- sortDFByGene(haplo.db$kval.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    count.df <- sortDFByGene(haplo.db$count.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    allele_palette <- alleleHapPalette(geno.df$ALLELES)
    AlleleCol <- allele_palette$AlleleCol
    transper <- allele_palette$transper

    ########################################################################################################

    ### Prepare All panels

    ## Middle panel
    p = ggplot(geno.df, aes(x = GENE, fill = factor(ALLELES,levels=AlleleCol))) + theme_bw() +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            text = element_text(size = text_size), strip.background = element_blank(),
            strip.text = element_text(face = "bold"),axis.text = element_text(colour = "black"),
            panel.spacing = unit(0, "cm"),
            strip.switch.pad.grid = unit(0, "cm"),
            plot.margin = unit(c(0.25, 0, 0.2, 0), "cm")) + geom_bar(position = "fill",width = 0.9,na.rm = T) +
      coord_flip() + xlab("") + ylab("") + facet_grid(paste0(".~", "hapBy"),switch = 'x') +
      scale_fill_manual(values=alpha(names(AlleleCol),transper),name='Alleles')

    if(!plotYaxis){
      p=p+theme(axis.text.y=element_blank())
    }

    ## Right panel
    ## plot K values
    pk <- ggplot(kval.df, aes(x = GENE, fill = K_GROUPED)) + theme_bw() +
      theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            text = element_text(size = text_size), strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(0, "cm"),
            strip.switch.pad.grid = unit(0, "cm"),
            plot.margin = unit(c(0.25, 0, 0.2, 0), "cm")) + geom_bar(position = "fill",width = 0.7,na.rm = T)  +
      coord_flip() +xlab("")+ ylab("") + facet_grid(paste0(".~", "hapBy"),switch = 'x')


    ## Left panel
    p2 <- ggplot(count.df,aes(x=GENE,y=COUNT2,fill=factor(ALLELES,levels=AlleleCol))) + geom_bar(stat="identity",position = 'Dodge',width = 0.9,na.rm = T) +
      coord_flip() + background_grid(minor='none')+scale_fill_manual(values=alpha(names(AlleleCol),transper),name='ALLELES') +
      theme(legend.position="none",strip.text = element_text(face = "bold"),axis.text = element_text(colour = "black"),
            text = element_text(size = text_size),plot.margin = unit(c(0.25, 0, -0.05, 0), "cm"),panel.background = element_blank()) + scale_y_continuous(breaks = seq(-3, 3, by = 1),labels = c(3:0,1:3)) +
      ylab(expression('log'[10]*'(Count+1)')) + xlab('Gene')  + geom_hline(yintercept=c(0), linetype="dotted")

    ########################################################################################################

    ### Plot All panels

    if(html_output){

      ## Prepare panels for html plot

      p = p + do.call(theme, list(...)) + theme( axis.title.x=element_blank())
      p.l <- ggplotly(p,height = 1000,width = 700) %>% plotly::layout(showlegend=FALSE)

      pk = pk + do.call(theme, list(...)) + theme( axis.title=element_blank())
      pk <- pk +   scale_fill_brewer(name='log<sub>10</sub>(lK)',drop = FALSE)
      pk.l <- ggplotly(pk,height = 1000,width = 700) %>% plotly::layout(showlegend=TRUE)
      pk.l$x$layout$annotations[[1]]$text = p.l$x$layout$annotations[[1]]$text
      pk.l$x$layout$annotations[[2]]$text = p.l$x$layout$annotations[[2]]$text

      p2 <- p2 + ylab('log<sub>10</sub>(Count+1)')
      p2.l <- ggplotly(p2,height = 1300,width = 1000) %>% plotly::layout(margin=list(b=50),
                                                                 yaxis = list(title = paste0(c(rep("&nbsp;", 3),
                                                                                               "Gene",
                                                                                               rep("&nbsp;", 3),
                                                                                               rep("\n&nbsp;", 1)),
                                                                                             collapse = "")),showlegend=TRUE)

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
          if(!is.na(NA)){
            count <- strsplit(strsplit(label,'<br />Count: ')[[1]][2],'<')[[1]][1]
            if(count=='NA') next
            count <- as.numeric(count)
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
      p.l.c$x$layout$annotations[[6]]$text = "log<sub>10</sub>(lK)"
      p.l.c$x$layout$annotations[[6]]$xanchor="center"
      p.l.c$x$layout$annotations[[6]]$y=0.99-0.0233*(length(AlleleCol)+2.4)#0.52
      p.l.c$x$layout$annotations[[6]]$x=0.98


      p.l.c$x$layout$annotations[[3]] <- p.l.c$x$layout$annotations[[6]]
      p.l.c$x$layout$annotations[[3]]$text = "Alleles"
      p.l.c$x$layout$annotations[[3]]$y=0.99
      p.l.c$x$layout$annotations[[3]]$x=0.98
      p.l.c$x$layout$annotations[[3]]$legendTitle=FALSE

      for(i in (length(AlleleCol)-2+1):(length(p.l.c$x$data)-(length(pk.l$x$data)+4))){
        p.l.c$x$data[[i]]$showlegend <- FALSE
      }

      plot_list[[sample_name]] <- p.l.c

    }else{
      p.legend <- get_legend(p)
      p = p +theme(legend.position='none')

      pk = pk + scale_fill_brewer(name=expression('log'[10]*'(lK)'),drop = FALSE)
      pk.legend <- get_legend(pk)
      pk = pk +theme(legend.position='none')

      p = p + do.call(theme, list(...)) + theme( axis.title.x=element_blank())
      pk = pk + do.call(theme, list(...)) + theme( axis.title=element_blank())

      p.legends  <-  plot_grid(pk.legend,p.legend, ncol=1, rel_heights = c(0.5,0.5),align = 'hv')

      p1 <- plot_grid(p2,p,pk,nrow=1,rel_widths = c(0.35,0.15,0.05),align = 'hv',axis='b')

      p <- plot_grid(p1,p.legends,ncol=2,rel_widths = c(1,0.1))
      #p1 <- ggdraw(plot_grid(p,pk,p.legends, nrow=1, rel_widths=c(0.1,0.1 ,0.1)))

      #p <- plot_grid(p2 + theme(plot.margin = unit(c(0.5,0,0,0),"lines")),
      #               p1 + theme(plot.margin = unit(c(0,0,0.5,0),"lines")),
      #               nrow=1,rel_widths = c(1,2))
      # now add the title
      #title <- ggdraw() + draw_label(sample_name, fontface='bold')

      plot_list[[sample_name]] <- p#plot_grid(title, p, ncol=1, rel_heights=c(0.05, 1))
      }

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
      for(sample_name in names(plot_list)){

        title <- ggdraw() + draw_label(sample_name, fontface='bold')
        plot(plot_grid(title, plot_list[[sample_name]], ncol=1, rel_heights=c(0.05, 1)))}

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
#' @param    kThreshDel           The minimum lK (log10 of the Bayes factor) used in \code{createFullHaplotype} to call a deletion. Indicates the color for strong deletion. Defualt is 3.
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
deletionHeatmap <- function(hap_table,html_output=FALSE,kTreshDel=3,chain=c('IGH','IGK','IGL'), ...){

  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(!("SUBJECT" %in% names(hap_table))){hap_table$SUBJECT <- rep('S1',nrow(hap_table))}

  hapBy_cols = names(hap_table)[grep(chain,names(hap_table))]

  genes_hap <- unique(substr(hap_table$GENE,4,4))

  GENE.loc.tmp <- GENE.loc[[chain]][grep(paste0(genes_hap,collapse = '|'),GENE.loc[[chain]])]

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

  hap_table.del.heatmap$DEL <- ifelse(hap_table.del.heatmap$ALLELE == 'Del'&hap_table.del.heatmap$K>=kTreshDel,3,0 )
  hap_table.del.heatmap$DEL[(! hap_table.del.heatmap$ALLELE  %in% c('Del','Unk','NA')) & hap_table.del.heatmap$K<3] <- 1
  hap_table.del.heatmap$DEL[hap_table.del.heatmap$ALLELE == 'Del'&hap_table.del.heatmap$K<kTreshDel] <- 2
  hap_table.del.heatmap$DEL[hap_table.del.heatmap$ALLELE == 'NA'] <- 4



  #manual reshape
  hap_table.del.heatmap.02 <- hap_table.del.heatmap[hap_table.del.heatmap$HapBy==ALLELE_01_num,]
  hap_table.del.heatmap.03 <- hap_table.del.heatmap[hap_table.del.heatmap$HapBy==ALLELE_02_num,]
  hap_table.del.heatmap.02$GENE2 <- gsub('IG[H|K|L]','',hap_table.del.heatmap.02$GENE)
  hap_table.del.heatmap.03$GENE2 <- gsub('IG[H|K|L]','',hap_table.del.heatmap.03$GENE)

  hap_table.del.heatmap.02$SUBJECT <- factor(x = hap_table.del.heatmap.02$SUBJECT,
                                             levels = unique(hap_table.del.heatmap.02$SUBJECT))

  hap_table.del.heatmap.02$GENE2 <- factor(x = hap_table.del.heatmap.02$GENE2,
                                           levels = gsub('IG[H|K|L]','',GENE.loc.tmp))

  hap_table.del.heatmap.03$SUBJECT <- factor(x = hap_table.del.heatmap.03$SUBJECT,
                                             levels = unique(hap_table.del.heatmap.03$SUBJECT))

  hap_table.del.heatmap.03$GENE2 <- factor(x = hap_table.del.heatmap.03$GENE2,
                                           levels = gsub('IG[H|K|L]','',GENE.loc.tmp))
  heatmap.df <- rbind(hap_table.del.heatmap.02,hap_table.del.heatmap.03)

  ALLELE_01_col = gsub('_','*',gsub('IG[H|K|L]','',ALLELE_01_col))
  ALLELE_02_col = gsub('_','*',gsub('IG[H|K|L]','',ALLELE_02_col))
  heatmap.df$DEL <- factor(heatmap.df$DEL,levels=c(0:4))
  heatmap.df$HapBy <- ifelse(heatmap.df$HapBy==ALLELE_01_num,ALLELE_01_col,ALLELE_02_col)
  heatmap.plot <- ggplot(data = heatmap.df, aes(x = GENE2, y = SUBJECT)) + theme_bw()+
    geom_tile(aes(fill = DEL)) + facet_wrap(~HapBy,nrow=2)+ scale_x_discrete(drop=FALSE)+
    scale_fill_manual(name='lK',labels=c('No deletion (lK>=3)','No deletion (lK<3)',paste0('Deletion (lK<',kTreshDel,')'),paste0('Deletion (lK>=',kTreshDel,')'),'NA'),values = c('white','lightgrey','lightblue','blue','gray40'),drop = FALSE)  + ylab('Subject') + xlab('Gene') +
    theme(strip.text = element_text(size=18),axis.title = element_text(size=18),axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90,vjust = 0.5 ,hjust = 1),plot.margin = margin(b = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.direction = "horizontal",
          legend.justification="center" ,
          legend.box.just = "bottom",
          legend.text=element_text(size=16),
          legend.key = element_rect(fill = "white", colour = "black"))


  del.df.heatmap <- heatmap.df
  del.df.heatmap <- del.df.heatmap %>% filter(ALLELE=='Del')

  del.df.heatmap.cnt <- del.df.heatmap %>% group_by(SUBJECT,GENE2) %>% mutate(n=n()) %>% dplyr::slice(1)
  del.df.heatmap.cnt$HapBy[del.df.heatmap.cnt$n==2] <- "Both"
  del.df.heatmap.cnt <- del.df.heatmap.cnt %>% group_by(GENE2,HapBy) %>% count_()

  del.df.heatmap.cnt$HapBy <- factor(del.df.heatmap.cnt$HapBy,levels=c(ALLELE_01_col,ALLELE_02_col,'Both'))
  pdel <- ggplot(del.df.heatmap.cnt,aes(x=GENE2,y=nn,fill=(HapBy))) + theme_bw() +
    geom_bar(stat='identity',position = 'stack',na.rm = T)  +
    theme(strip.background = element_blank(),
          axis.text = element_text(size=14),
          axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),plot.margin = margin(0,8,0,7,"pt"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    ylab('Number of individuals\nwith a deletion')+
    scale_fill_manual(name = "Chromosome",values = c('darksalmon','deepskyblue','darkolivegreen3','grey50'))  +
    scale_x_discrete(drop=FALSE) + xlab('')

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



  comb <- plot_grid(pdel, heatmap.plot,ncol=1,rel_heights=c(0.15, 0.3),align = "hv",axis='b')
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


  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  hapBy_alleles <-  gsub('_','*',names(hap_table)[grep(chain,names(hap_table))])
  hapBy_cols <- gsub('IG[H|K|L]','',hapBy_alleles)

  haplo_db <- c()
  for(sample_name in unique(hap_table$SUBJECT)){
    geno.df <- parseHapTab(hap_table[hap_table$SUBJECT==sample_name,], chain = chain, df_ToReturn = 'geno.df')
    geno.df <- sortDFByGene(geno.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    geno.df$SUBJECT <- sample_name
    haplo_db <- rbind(haplo_db, geno.df)
  }

  allele_palette <- alleleHapPalette(haplo_db$ALLELES)

  heatmap.df <- haplo_db %>% group_by(SUBJECT,hapBy,GENE) %>% mutate(n=n())
  heatmap.df$freq <-ifelse(heatmap.df$n==2,0.5,1)
  heatmap.df$GENE <- factor(heatmap.df$GENE, levels =  gsub('IG[H|K|L]','',GENE.loc[[chain]]))
  heatmap.df$title <- ifelse(heatmap.df$hapBy==hapBy_cols[1],hapBy_cols[1],hapBy_cols[2])

  col <- ifelse(c(1:(nrow(unique(heatmap.df[heatmap.df$hapBy==hapBy_cols[1],'SUBJECT']))*4))%%3==0,'black','white')

  p <- ggplot(heatmap.df[heatmap.df$hapBy==hapBy_cols[1],], (aes(x = GENE, y = freq, fill = factor(ALLELES,levels=allele_palette$AlleleCol)))) +
    geom_col(position = "fill", width = 0.95,na.rm = T) +
    scale_fill_manual(values=alpha(names(allele_palette$AlleleCol),allele_palette$transper),name='Alleles',drop=F)+
    facet_grid(SUBJECT ~ title, as.table = FALSE, switch = "y") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +  theme(axis.text.x = element_text(angle = 90,vjust = 0.5 ,hjust = 1, size=10,colour = 'black'),
                                                 strip.text.x = element_text(size=14),
                                                 strip.text.y = element_text(angle = 180, size=10),
                                                 panel.grid = element_blank(),
                                                 strip.placement = "outside",
                                                 axis.ticks.y = element_line(colour = col),
                                                 axis.text.y = element_blank(),
                                                 strip.background.y = element_blank(),
                                                 panel.spacing.y = unit(0.9, "pt")) + labs(y='',x='')

  col <- ifelse(c(1:(nrow(unique(heatmap.df[heatmap.df$hapBy==hapBy_cols[2],'SUBJECT']))*4))%%3==0,'black','white')

  p1 <- ggplot(heatmap.df[heatmap.df$hapBy==hapBy_cols[2],], (aes(x = GENE, y = n, fill = factor(ALLELES,levels=allele_palette$AlleleCol)))) +
    geom_col(position = "fill", width = 0.95,na.rm = T) +
    scale_fill_manual(values=alpha(names(allele_palette$AlleleCol),allele_palette$transper),name='Alleles',drop=F)+
    facet_grid(SUBJECT ~ title, as.table = FALSE, switch = "y") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +  theme(axis.text.x = element_text(angle = 90,vjust = 0.5 ,hjust = 1, size=10,colour = 'black'),
                                                 strip.text.x = element_text(size=14),
                                                 strip.text.y = element_text(angle = 180, size=10),
                                                 panel.grid = element_blank(),
                                                 strip.placement = "outside",
                                                 axis.ticks.y = element_line(colour = col),
                                                 axis.text.y = element_blank(),
                                                 strip.background.y = element_blank(),
                                                 panel.spacing.y = unit(0.9, "pt")) + labs(y='')

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
  GENE.loc.tmp <- GENE.loc[[chain]][GENE.loc[[chain]] %in% GENE.usage.df$GENE]

  GENE.usage.df$GENE2 <- factor(gsub(chain,'',GENE.usage.df$GENE), levels=gsub(chain,'',GENE.loc.tmp))

  colvec <- ifelse(GENE.loc.tmp%in%genes.low.cer, "red", ifelse(GENE.loc.tmp%in%genes.dup, "purple", "black"))

  ### gene usage with deletions in population according to binom test
  p.del <- ggplot(GENE.usage.df %>% filter(DELETION!='Non reliable'),aes(x=GENE2,y=FRAC)) + geom_boxplot(outlier.colour=NA) +
    geom_jitter(aes(x=GENE2,color=DELETION),width = 0.25,size=0.5)+ theme(axis.text.y = element_text(size=16),
                                                                     axis.title = element_text(size=16),
                                                                     axis.text.x = element_text(size=14,angle = 90, hjust = 1,vjust=0.5,color=colvec),
                                                                     legend.text = element_text(size=16),
                                                                     legend.position = 'none')+
    ylab('Fraction') + xlab('')  + scale_color_manual(name='',labels=c('Deletion','No Deletion','NA'),values = c('blue','black','grey40'),drop=T)+
    guides(color = guide_legend(override.aes = list(size=5)))

  ### heat map of deletions in population according to binom test
  GENE.usage.df$DELETION <- factor(GENE.usage.df$DELETION,levels=levels(GENE.usage.df$DELETION))
  if(length(levels(GENE.usage.df$DELETION))<4) lab = c('Deletion','No Deletion','NA') else lab = c('Deletion','No Deletion','NA','Non reliable')
  if(length(levels(GENE.usage.df$DELETION))<4) col_val = c('#0000ff','#ffffff','#808080') else col_val = c('#0000ff','#ffffff','#808080','#ffefd5')
  heatmap.plot <-   ggplot(data = GENE.usage.df, aes(x = GENE2, y = SUBJECT)) +
    geom_tile(aes(fill = DELETION)) +
    scale_fill_manual(name='',labels=lab,values = col_val,drop=F)+
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


########################################################################################################
#' Hierarchical clustering of haplotypes graphical output
#'
#' \code{hapDendo} creates a graphical output of an hierarchical clustering based on the Jaccardian distance between multiple individuals haplotypes.
#'
#' @details A \code{data.frame} created by \code{createFullHaplotype}.
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    gene_sort            If by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names. Defualt is TRUE.
#' @param    mark_low_lk          if TRUE, a texture is add for low lK values. Defualt is TRUE.
#' @param    lk_cutoff            the lK cutoff value to be considerd low for texture layer. Defualt is lK<1.
#'
#' @return   a single chromosome deletion visualization.
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db, HVGERM, HDGERM)
#'
#' # Infering haplotype
#' hap_df = createFullHaplotype(samples_db,toHap_col=c("V_CALL","D_CALL"),hapBy_col="J_CALL",hapBy="IGHJ6",toHap_GERM=c(HVGERM,HDGERM));
#' hapDendo(hap_df)
#'
#' @export
hapDendo <- function(hap_table, chain = c('IGH','IGK','IGL'), gene_sort = c("name", "position"), removeIGH=TRUE, mark_low_lk=TRUE, lk_cutoff=1,...){

  if(missing(chain)) {
    chain='IGH'
  }
  chain <- match.arg(chain)

  if(missing(gene_sort)) {
    gene_sort='position'
  }
  gene_sort <- match.arg(gene_sort)

  hapBy_cols = names(hap_table)[grep(chain,names(hap_table))]
  samples <- unique(hap_table$SUBJECT)

  # creating the distance matrix for clustering
  mat <- matrix(NA,nrow = length(samples),ncol = length(samples))
  for(i in 2:length(samples)){
    for(j in 1:(i-1)){

      hap_merge <- merge(hap_table[hap_table$SUBJECT==samples[i], c("GENE",hapBy_cols[1],hapBy_cols[2])],
                         hap_table[hap_table$SUBJECT==samples[j], c("GENE",hapBy_cols[1],hapBy_cols[2])],
                         by='GENE')

      mat[i,j] <-  calcJacc(vec1A = hap_merge[,2],vec1B = hap_merge[,3],
                            vec2A = hap_merge[,4],vec2B = hap_merge[,5], method = 'geneByGene')

    }
  }
  colnames(mat) <- samples
  rownames(mat) <- samples

  # finding the hierarchical clustering
  fit <- hclust(as.dist(mat), method="ward.D")
  samples <- samples[fit$order]

  # Preparing the hclust data for plotting
  dend <- as.dendrogram(hclust(as.dist(mat), method="ward.D"))
  dend_data <- ggdendro::dendro_data(dend)
  segment_data <- with(ggdendro::segment(dend_data), data.frame(x = y, y = x, xend = yend, yend = xend))

  # Using the dendrogram label data to position the samples labels
  samples_pos_table <- with(dend_data$labels,data.frame(y_center = x, gene = as.character(label), height = 1))
  samples_axis_limits <- with(samples_pos_table, c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1)

  plt_dendr <- ggplot(segment_data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(breaks = samples_pos_table$y_center,
                       labels = samples_pos_table$gene,
                       limits = samples_axis_limits,
                       expand = c(0, 0)) +
    labs(x = "Jaccardian distance", y = "", colour = "", size = "") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.ticks.y =  element_blank(),
          panel.grid = element_blank(),panel.border = element_blank(),
          axis.text = element_text(size=10,colour = 'black'),
          axis.title.x = element_text(size=14,colour = 'black'))

  # Creating the haploype db for ploting
  haplo_db_clust <- c()
  for(sample_name in samples){
    geno.df <- parseHapTab(hap_table[hap_table$SUBJECT==sample_name,], chain = chain, df_ToReturn = 'geno.df')
    kval.df <- parseHapTab(hap_table[hap_table$SUBJECT==sample_name,], chain = chain, df_ToReturn = 'kval.df')
    geno.df <- sortDFByGene(geno.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    kval.df <- sortDFByGene(kval.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    geno.df$K <- apply(geno.df[,c('GENE','hapBy')],1,function(x){kval.df$K[kval.df$GENE==x[[1]]&kval.df$hapBy==x[[2]]]})
    geno.df$SUBJECT <- sample_name
    haplo_db_clust <- rbind(haplo_db_clust, geno.df)
  }
  allele_palette <- alleleHapPalette(haplo_db_clust$ALLELES)

  # Formating the data to fit heatmap plot
  haplo_db_clust <- haplo_db_clust %>% group_by(SUBJECT,hapBy,GENE) %>% mutate(n=n())
  haplo_db_clust$freq <-ifelse(haplo_db_clust$n==2,0.5,1)
  haplo_db_clust$GENE <- factor(haplo_db_clust$GENE, levels =  gsub('IG[H|K|L]','',GENE.loc[[chain]]))
  haplo_db_clust$grouper_x <- 'Gene'
  haplo_db_clust$grouper_y <- sapply(1:nrow(haplo_db_clust),function(i) paste0(haplo_db_clust$SUBJECT[i],' ',haplo_db_clust$hapBy[i]))

  haplo_db_clust_texture <- c()
  loc <- 1:length(levels(droplevels(haplo_db_clust$GENE)))
  names(loc) <- levels(droplevels(haplo_db_clust$GENE))
  for(i in 1:nrow(haplo_db_clust)){
    if(haplo_db_clust$K[i]<lk_cutoff && !haplo_db_clust$ALLELES[i] %in% c('Unk','Del','NR')){
      tmp_point <- haplo_db_clust[i,] %>% slice(rep(1,each=ifelse(length(samples)<4,15,8))) %>% mutate(points=seq(0,0.9,length.out = ifelse(length(samples)<4,15,8)),yend=seq(0,0.9,length.out = ifelse(length(samples)<4,15,8))+0.1,GENE2=loc[as.character(GENE)])
      haplo_db_clust_texture <- rbind(haplo_db_clust_texture,tmp_point)
    }
  }
  haplo_db_clust_texture <- haplo_db_clust_texture[!duplicated(haplo_db_clust_texture[,c(1,2,4,10)]),]
  allele_cols <- gsub('IG[H|K|L]','',gsub('_','*',hapBy_cols))
  allele_palette <- alleleHapPalette(haplo_db_clust$ALLELES)

  # Adding white space to plot
  heatmap.df <- c()
  samples_order <- c()
  samples_label <- c()
  for(i in 1:length(samples)){
    samp = samples[i]
    sub <- haplo_db_clust[haplo_db_clust$SUBJECT==samp,]
    sub2 <- sub[sub$hapBy==allele_cols[1],]
    sub2$freq <- 0
    sub2$grouper_y <- paste0(sub2$SUBJECT[1],' NA')
    if(i != length(samples)){
      sub <- rbind(sub,sub2)
      heatmap.df <- rbind(heatmap.df,sub)
      samples_order <- c(samples_order,unique(sub$grouper_y))
      tmp_l <- c(unique(sub$hapBy),'')
      names(tmp_l) <- unique(sub$grouper_y)
      samples_label <- c(samples_label,tmp_l)
    }else{
      heatmap.df <- rbind(heatmap.df,sub)
      samples_order <- c(samples_order,unique(sub$grouper_y))
      tmp_l <- unique(sub$hapBy)
      names(tmp_l) <- unique(sub$grouper_y)
      samples_label <- c(samples_label,tmp_l)}
  }

  heatmap.df$grouper_y <- factor(heatmap.df$grouper_y,levels=samples_order)
  haplo_db_clust_texture$grouper_y <- factor(haplo_db_clust_texture$grouper_y,levels=samples_order)
  hap_plot <- ggplot(heatmap.df, (aes(x = GENE, y = freq, fill = factor(ALLELES,levels=allele_palette$AlleleCol)))) +
    geom_col(position = "fill", width = 0.95,na.rm = T) +
    scale_fill_manual(values=alpha(names(allele_palette$AlleleCol),
                                   allele_palette$transper),
                      name='Alleles',drop=F)+
    facet_grid(grouper_y~grouper_x  , as.table = FALSE, switch = "y",labeller=labeller(grouper_y=samples_label)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +  theme(axis.text.x = element_text(
      angle = 90,vjust = 0.5 ,hjust = 1, size=10,colour = 'black'),
      strip.text.x = element_blank(),
      strip.text.y = element_text(angle = 180, size=10),
      panel.grid = element_blank(),
      strip.placement = "outside",
      axis.ticks.y = element_line(colour = 'white'),
      axis.line.y.left = element_blank(),
      axis.text.y = element_blank(),
      strip.background.y = element_blank(),
      strip.background.x = element_blank(),
      panel.spacing.y = unit(0.9, "pt"),
      legend.position="bottom",
      axis.title.x = element_text(size=14,colour = 'black'),
      legend.justification="center") + labs(y='',x='Gene') +
    guides(fill = guide_legend(nrow = round(length(allele_palette$AlleleCol)/7),order = 1))



  if(mark_low_lk){
    # Get Allele legend
    gt1 = ggplotGrob(hap_plot)

    hap_plot <- hap_plot + geom_segment(data=haplo_db_clust_texture,aes(x=GENE2-.49, xend=GENE2+.49, y=points, yend=yend,color='<1'))+
    scale_color_manual(values = c('white'),name='lK') +
    guides(color=guide_legend(override.aes = list(size = 1),order = 2),fill="none") +
    theme(legend.justification="center",legend.key = element_rect(fill = 'gray'))


    # Get lK legend
    gt2 = ggplotGrob(hap_plot)

    leg1 = gtable::gtable_filter(gt1, "guide-box")
    leg2 = gtable::gtable_filter(gt2, "guide-box")
    # Combine the legends
    leg <- cbind(leg1[["grobs"]][[1]],  leg2[["grobs"]][[1]], size = "first")
    # Insert legend into g1 (or g2)
    gt1$grobs[gt1$layout$name == "guide-box"][[1]] <- leg
    gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[3,c('t','b')] <- gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[1,c('t','b')]
    gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[3,c('l','r')] <- gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[1,c('l','r')] + 2
    legend <- get_legend(gt1)
    } else legend <- get_legend(hap_plot)

  dend_hap <- plot_grid(hap_plot+theme(legend.position = 'none'),plt_dendr,nrow=1,align = 'h',axis = 'b',rel_widths = c(1,0.25))
  plot(plot_grid(dend_hap,legend,ncol=1,rel_heights = c(1,0.1)))
}
