# Generate sysdata

#### Default Gene location ####

## Get the tables from imgt
imgt_igk = "http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=GeneOrder&species=human&group=IGK"
imgt_igk  = XML::readHTMLTable(imgt_igk, header=T, which=1,stringsAsFactors=F, skip.rows = 1)
cat(paste0('"',rev(imgt_igk$`IMGT gene name`),'"',collapse = ','))

imgt_igL = "http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=GeneOrder&species=human&group=IGL"
imgt_igL  = XML::readHTMLTable(imgt_igL, header=T, which=1,stringsAsFactors=F, skip.rows = 1)
cat(paste0('"',rev(imgt_igL$`IMGT gene name`),'"',collapse = ','))

imgt_igH = "http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=GeneOrder&species=human&group=IGH"
imgt_igH  = XML::readHTMLTable(imgt_igH, header=T, which=1,stringsAsFactors=F, skip.rows = 1)
cat(paste0('"',rev(imgt_igH$`IMGT gene name`),'"',collapse = ','))

GENE.loc <- list(

  IGH = c("IGHV7-81", "IGHV5-78", "IGHV3-74", "IGHV3-73", "IGHV3-72", "IGHV3-71","IGHV2-70", "IGHV1-69D", "IGHV1-f","IGHV2-70D", "IGHV1-69", "IGHV1-68","IGHV1-69-2", "IGHV3-69-1",
          "IGHV3-66","IGHV3-64", "IGHV3-63", "IGHV3-62", "IGHV4-61" , "IGHV4-59","IGHV1-58", "IGHV4-55", "IGHV3-54", "IGHV3-53", "IGHV3-52",
          "IGHV5-51", "IGHV3-49", "IGHV3-48", "IGHV3-47", "IGHV1-46", "IGHV1-45", "IGHV3-43" , "IGHV3-43D", "IGHV7-40", "IGHV4-39", "IGHV1-38-4", "IGHV3-38-3", "IGHV4-38-2",
          "IGHV3-38", "IGHV3-35", "IGHV7-34-1", "IGHV4-34", "IGHV3-33-2", "IGHV3-33", "IGHV3-32",
          "IGHV4-31", "IGHV3-30-52" , "IGHV3-30-5", "IGHV3-30-42", "IGHV4-30-4", "IGHV3-30-33", "IGHV3-30-3", "IGHV3-30-22", "IGHV4-30-2", "IGHV4-30-1", "IGHV3-30-2", "IGHV3-30", "IGHV3-29", "IGHV4-28" ,
          "IGHV2-26", "IGHV3-25" ,"IGHV1-24", "IGHV3-23", "IGHV3-23D", "IGHV3-22", "IGHV3-21","IGHV3-20", "IGHV3-19","IGHV1-18", "IGHV3-16",
          "IGHV3-15","IGHV3-13", "IGHV3-11" , "IGHV2-10", "IGHV3-9","IGHV5-10-1" , "IGHV1-8","IGHV3-64D", "IGHV3-7",
          "IGHV2-5","IGHV7-4-1","IGHV4-4", "IGHV1-3", "IGHV1-2",  "IGHV6-1",
          "IGHD1-1" ,  "IGHD2-2" , "IGHD3-3", "IGHD6-6" , "IGHD1-7" ,  "IGHD2-8" ,
          "IGHD3-9",  "IGHD3-10"  ,  "IGHD4-11", "IGHD5-12" , "IGHD6-13", "IGHD1-14",
          "IGHD2-15","IGHD3-16",
          "IGHD4-17","IGHD5-18"   , "IGHD6-19","IGHD1-20"  ,"IGHD2-21"  ,"IGHD3-22" ,
          "IGHD4-23" ,  "IGHD5-24" ,  "IGHD6-25"  , "IGHD1-26","IGHD7-27",
          "IGHJ1", "IGHJ2","IGHJ3","IGHJ4","IGHJ5","IGHJ6"),

  # IGK =  c("IGKJ5","IGKJ4","IGKJ3","IGKJ2","IGKJ1","IGKV4-1","IGKV5-2","IGKV1-5","IGKV1-6","IGKV1-8","IGKV1-9","IGKV3-11",
  #          "IGKV1E-12","IGKV1E-13","V3-15","IGKV1-16","IGKV1-17","IGKV3-20","IGKV6E-21","IGKV2-24","IGKV1-27","IGKV2E-28","IGKV2-29","IGKV2-30",
  #          "IGKV1E-33","IGKV1E-39","V2E-40","IGKV2D-30","IGKV2D-29",
  #          "IGKV2D-26","IGKV3D-20","V1D-17","IGKV1D-16","IGKV3D-15","IGKV3D-11","IGKV1D-43","IGKV1D-8","IGKV3D-7"),

  IGK =  c("IGKJ5","IGKJ4","IGKJ3","IGKJ2","IGKJ1","IGKV4-1","IGKV5-2","IGKV7-3","IGKV2-4","IGKV1-5","IGKV1-6","IGKV3-7",
           "IGKV1-8","IGKV1-9","IGKV2-10","IGKV3-11","IGKV1-12","IGKV1-13","IGKV2-14","IGKV3-15","IGKV1-16","IGKV1-17","IGKV2-18",
           "IGKV2-19","IGKV3-20","IGKV6-21","IGKV1-22","IGKV2-23","IGKV2-24","IGKV3-25","IGKV2-26","IGKV1-27","IGKV2-28","IGKV2-29",
           "IGKV2-30","IGKV3-31","IGKV1-32","IGKV1-33","IGKV3-34","IGKV1-35","IGKV2-36","IGKV1-37","IGKV2-38","IGKV1-39","IGKV2-40",
           "IGKV2D-40","IGKV1D-39","IGKV2D-38","IGKV1D-37","IGKV2D-36","IGKV1D-35","IGKV3D-34","IGKV1D-33","IGKV1D-32","IGKV3D-31",
           "IGKV2D-30","IGKV2D-29","IGKV2D-28","IGKV1D-27","IGKV2D-26","IGKV3D-25","IGKV2D-24","IGKV2D-23","IGKV1D-22","IGKV6D-21",
           "IGKV3D-20","IGKV2D-19","IGKV2D-18","IGKV6D-41","IGKV1D-17","IGKV1D-16","IGKV3D-15","IGKV2D-14","IGKV1D-13","IGKV1D-12",
           "IGKV3D-11","IGKV2D-10","IGKV1D-42","IGKV1D-43","IGKV1D-8","IGKV3D-7","IGKV1-NL1"),

  IGL =  c("IGLV4-69","IGLV8-61","IGLV4-60","IGLV6-57","IGLV10-54","IGLV5-52","IGLV1-51","IGLV9-49",
           "IGLV1-47","IGLV7-46","IGLV5-45","IGLV1-44","IGLV7-43","IGLV1-40","IGLV5-37","IGLV1-36",
           "IGLV3-27","IGLV3-25","IGLV2-23","IGLV3-22","IGLV3-21","IGLV3-19","IGLV2-18","IGLV3-16","IGLV2-14",
           "IGLV3-12","IGLV2-11","IGLV3-10","IGLV3-9","IGLV2-8","IGLV4-3","IGLV3-1","IGLJ1","IGLJ2","IGLJ3","IGLJ6","IGLJ7")

)

#### Known Pseudo Gene####

library(rvest)
imgt = "http://www.imgt.org/genedb/GENElect?query=4.2+IG&species=Homo+sapiens"
pseudo <- html_nodes(webpage, "table") %>% html_table()
pseudo_IGK <- pseudo[[1]] %>% filter(`IMGT gene functionality`%in% c("ORF", "P"), grepl('IGK', `IMGT/GENE-DB`))
cat(paste0('"',rev(pseudo_IGK$`IMGT/GENE-DB`),'"',collapse = ','))

pseudo_IGL <- pseudo[[1]] %>% filter(`IMGT gene functionality`%in% c("ORF", "P"), grepl('IGL', `IMGT/GENE-DB`))
cat(paste0('"',rev(pseudo_IGL$`IMGT/GENE-DB`),'"',collapse = ','))


PSEUDO <- list(

  IGH = c("IGHV2-10","IGHV3-52","IGHV3-47","IGHV3-71","IGHV3-22","IGHV4-55","IGHV1-68","IGHV2-10","IGHV5-78","IGHV3-32","IGHV3-33-2","IGHV3-38-3",
          "IGHV3-25","IGHV3-19","IGHV7-40","IGHV3-63","IGHV3-62","IGHV3-29","IGHV3-54","IGHV1-38-4","IGHV7-34-1","IGHV1-38-4","IGHV3-30-2","IGHV3-69-1","IGHV3-30-22",
          "IGHV1-f","IGHV3-30-33","IGHV3-38","IGHV7-81","IGHV3-35","IGHV3-16","IGHV3-30-52","IGHV1-69D", "IGHD1-14", "IGHV3-30-42"),

  # IGK = c("IGKV1-22","IGKV1-35","IGKV1-39","IGKV1D-22","IGHV3-22","IGKV1D-27","IGKV1D-32","IGKV1D-35","IGKV2-10","IGKV2-14","IGKV2-18","IGKV2-19",
  #         "IGKV2-23","IGKV2-26","IGKV2-36","IGKV2-38","IGKV2-4","IGKV2D-10","IGKV2D-14","IGKV2D-18","IGKV2D-19","IGKV2D-23","IGKV2D-24","IGKV1D-42","IGKV2D-36","IGKV2D-38","IGKV3-25",
  #         "IGKV3-31","IGKV3-34","IGKV3-7","IGKV3D-25","IGKV3D-31","IGKV3D-34","IGKV6D-41","IGKV7-3","IGKV1D-42","IGKV1D-37","IGKV1E-37"),

  IGK = c("IGKV7-3","IGKV6D-41","IGKV3D-34","IGKV3D-31","IGKV3D-25","IGKV3/OR22-2","IGKV3/OR2-5","IGKV3/OR2-268","IGKV3-34","IGKV3-31","IGKV3-25",
          "IGKV2D-38","IGKV2D-36","IGKV2D-24","IGKV2D-23","IGKV2D-19","IGKV2D-18","IGKV2D-14","IGKV2D-10","IGKV2/OR22-4","IGKV2/OR22-3","IGKV2/OR2-8",
          "IGKV2/OR2-7D","IGKV2/OR2-7","IGKV2/OR2-4","IGKV2/OR2-2","IGKV2/OR2-10","IGKV2/OR2-1","IGKV2-4","IGKV2-38","IGKV2-36","IGKV2-26","IGKV2-23",
          "IGKV2-19","IGKV2-18","IGKV2-14","IGKV2-10","IGKV1D-42","IGKV1D-37","IGKV1D-35","IGKV1D-32","IGKV1D-27","IGKV1D-22","IGKV1/ORY-1",
          "IGKV1/OR9-2","IGKV1/OR9-1","IGKV1/OR22-5","IGKV1/OR22-1","IGKV1/OR2-9","IGKV1/OR2-6","IGKV1/OR2-3","IGKV1/OR2-2","IGKV1/OR2-118",
          "IGKV1/OR2-11","IGKV1/OR2-108","IGKV1/OR2-1","IGKV1/OR2-0","IGKV1/OR15-118","IGKV1/OR10-1","IGKV1/OR1-1","IGKV1/OR-4","IGKV1/OR-3",
          "IGKV1/OR-2","IGKV1-37","IGKV1-35","IGKV1-32","IGKV1-22"),

  IGL = c("IGLV7-35","IGLV3-7","IGLV3-6","IGLV3-4","IGLV3-32","IGLV3-31","IGLV3-30","IGLV3-29","IGLV3-26","IGLV3-24","IGLV3-2","IGLV3-17","IGLV3-15",
          "IGLV3-13","IGLV2-NL1","IGLV2-5","IGLV2-34","IGLV2-33","IGLV2-28","IGLV11-55","IGLV10-67","IGLV1-62","IGLV1-50","IGLV(VII)-41-1",
          "IGLV(VI)-25-1","IGLV(VI)-22-1","IGLV(V)-66","IGLV(V)-58","IGLV(IV)/OR22-2","IGLV(IV)/OR22-1","IGLV(IV)-66-1","IGLV(IV)-65","IGLV(IV)-64",
          "IGLV(IV)-59","IGLV(IV)-53","IGLV(I)-70","IGLV(I)-68","IGLV(I)-63","IGLV(I)-56","IGLV(I)-42","IGLV(I)-38","IGLV(I)-20","IGLL4","IGLL2",
          "IGLJ5","IGLJ4","IGLC5","IGLC4","IGLC/OR22-2","IGLC/OR22-1")
)

#### Default Haplotype Allele color scheme ####
ALLELE_PALETTE <- c(
  '01'='#f5bc6e','02'='#9d69f4','03'='#598200','04'='#e375fd','05'='#01a452',
  '06'='#ff6ed5','07'='#60de8a','08'='#e3006d','09'='#0edec9','10'='#c43e00',
  '11'='#0197f3','12'='#ffa629','13'='#7c82ff','14'='#8dd979','15'='#9a0f79',
  '16'='#008844','17'='#ff64aa','18'='#ff5d68','19'='#ff88ce','NA'='#00633c',
  'NA'='#355d06','NA'='#009f8b','NA'='#ff7d4e','NA'='#02a8ec','NA'='#7d4408',
  'NA'='#89c8ff','NA'='#ff897a','NA'='#58d8e0','NA'='#ff8793','NA'='#82d89e','NA'='#803a6a'
)

#### Default IGH Germlines ####

HVGERM <- toupper(unlist(seqinr::read.fasta('data-raw/HVGERM_4_12_18.fasta',as.string = T)))
HDGERM <- toupper(unlist(seqinr::read.fasta('data-raw/HDGERM_4_12_18.fasta',as.string = T)))
HJGERM <- toupper(unlist(seqinr::read.fasta('data-raw/HJGERM_4_12_18.fasta',as.string = T)))

save(HVGERM, file='data/HVGERM.rda')
save(HDGERM, file='data/HDGERM.rda')
save(HJGERM, file='data/HJGERM.rda')

#### Default Binom test cutoffs ####

Binom.test.gene.cutoff <- list(
  IGH = read.delim('data-raw/IGH_Binom_test_cutoff.tab',sep='\t',stringsAsFactors = FALSE)
)

#### Save to R/sysdata.rda ####

usethis::use_data(GENE.loc,
                   PSEUDO,
                   HVGERM,
                   HDGERM,
                   HJGERM,
                   Binom.test.gene.cutoff,
                   ALLELE_PALETTE,
                   internal=TRUE, overwrite=TRUE)

