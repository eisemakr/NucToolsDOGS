# This script takes the output of 2 DOGS2genes.R analysis detects which one belongs to mouse and which to men
# and compares which genes are under influence in both datasets by accessing either a file I downloaded manually
# or BiomaRt. Everything is set up for grch37 the code has to be modified if a different human genome would be used. 

#################################################################################

#Load packages, take inputs, check inputs, set wd, check wd

#################################################################################

if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.4.6.tar.gz",repos=NULL, method="libcurl")
pacman::p_load( biomaRt, dplyr )

args <- commandArgs(TRUE)
p2datahuman <- as.character(args[1])   # complete path to human DOGS2dirgene file
p2datamouse <- as.character(args[2])   # complete path to mouse DOGS2dirgene file
wd <- as.character(args[3])   # working directory. A Folder for the results will be created. 

#################################################################################
#for development
#################################################################################

#p2datahuman <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/Control/con1_con2/DOGS2dirgene1000.bed'
#p2datamouse <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/cutoff30/DOGS2dirgene1000.bed'
#wd <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/Control/dev'

#################################################################################

if(is.na(args[1])|is.na(args[2])|is.na(args[3])){
  print('please specify all 3 arguments!')
  quit()
}

setwd(wd)

if (getwd()!= wd) {
  print('working directory could not be set. Check pls thats essential!')
  quit()
}

#################################################################################
#read inputs and sort by Ensembl transcript ID 
#################################################################################

datahuman <- read.table(p2datahuman)
datamouse <- read.table(p2datamouse)

datahuman <- datahuman[order(datahuman$V4),]
datamouse <- datamouse[order(datamouse$V4),]


print('using grch37!!! Change my code if you need the latest human genome.')
print('when using grch38 for merging of data a getLDS function should be used.')
print('The code contains already some commented snippets for adjustion.')

#################################################################################
#Define Mart objects and accessing the corresponding data. 
#################################################################################

ensembl.human = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        host="grch37.ensembl.org",
                        path="/biomart/martservice",
                        dataset="hsapiens_gene_ensembl")

#ensembl.human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl.mouse= useMart("ensembl",dataset="mmusculus_gene_ensembl")

idhumanm <- getBM(attributes = c('mmusculus_homolog_associated_gene_name', 'mmusculus_homolog_chromosome',
                                 'mmusculus_homolog_chrom_start','mmusculus_homolog_chrom_end',
                                 'ensembl_transcript_id'), 
                 filters = 'ensembl_transcript_id',
                 values = datahuman$V4,
                 mart = ensembl.human)

idhuman <- getBM(attributes = c('chromosome_name', 'transcript_start','transcript_end',
                                'ensembl_transcript_id','hgnc_symbol'), 
                  filters = 'ensembl_transcript_id',
                  values = datahuman$V4,
                  mart = ensembl.human)
#Printing this dataframe should show the consistency of all downloaded data. 
#idtotalhuman <- cbind(datahuman[,c(1:4)],idhuman,idhumanm)

idmouseh <- getBM(attributes = c('hsapiens_homolog_associated_gene_name', 'hsapiens_homolog_chromosome',
                                 'hsapiens_homolog_chrom_start','hsapiens_homolog_chrom_end',
                                 'ensembl_transcript_id'), 
                 filters = 'ensembl_transcript_id',
                 values = datamouse$V4,
                 mart = ensembl.mouse)
idmouse <- getBM(attributes = c('chromosome_name', 'transcript_start','transcript_end',
                                'ensembl_transcript_id','mgi_symbol'), 
                 filters = 'ensembl_transcript_id',
                 values = datamouse$V4,
                 mart = ensembl.mouse)
#Printing this dataframe should show the consistency of all downloaded data. 
#idtotalmouse <- cbind(datamouse[,c(1:4)],idmouse,idmouseh)

#################################################################################
#Merging the data from both sides. 
#################################################################################

resmouse <- merge(idmouse,idhumanm, by.y = 'mmusculus_homolog_associated_gene_name', by.x = 'mgi_symbol')
resmouse <- resmouse[-which(resmouse$mgi_symbol==''),]
reshuman <- merge(idhuman,idmouseh, by.y = 'hsapiens_homolog_associated_gene_name', by.x = 'hgnc_symbol')
reshuman <- reshuman[-which(reshuman$hgnc_symbol==''),]

#################################################################################
#Printing the data to a txt file in the working directory. 
#################################################################################

sink('Geneoverlap.txt',append = T)
print(Sys.Date())
print(Sys.time())
print('inputs were:')
print(paste(as.character(args[1]), as.character(args[2]), as.character(args[3])))
print(resmouse)
print(reshuman)
sink()

quit()

#A Downloaded list of homologous mouse/human genes. It works also by extracting these values
anno <- read.table('X:/db05/Eike/mouse_human-txt.txt', fill = T, sep = '\t', header = T)

#This function snippet can later be used to make the join not by downlaoding 2 dataframes and then joining here but
#directly making a join db statement. Mistakes should then by on the site of ENSEMBL not us. 
#my.mouse.genes <- c('Frem2', 'Kmt2d', 'Scn8a', 'Abcg1', 'Acvr1b', 'Flnc', 'Lama2', 'Myh7b', 'Myo10', 'Ryr3')
#getLDS(attributes = c('mgi_symbol', 'chromosome_name'),
       #filters = 'mgi_symbol', values = my.mouse.genes, mart = ensembl.mouse,
       #attributesL = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position') , martL = ensembl.human)

#The old list of identifierst
#data <- read.table('X:/db05/Eike/genesmuethe_mirsi.txt', fill = T, sep = '\t')
#res <- vector(length = nrow(idmouse) )

#Extraction by list ... works fine
tt <- unique(anno[match(idmouse$mgi_symbol,anno$Symbol),1])
ttt <- unique(anno[match(idhuman$hgnc_symbol,anno$Symbol),1])
ttt[match(tt,ttt)]
anno[which(anno$HomoloGene.ID=='82993'),]

#This could acctually be deleted
# for (i in 1: nrow(idmouse)){
#   a<-0
#   print(i)
#   a<- grep(idmouse$mgi_symbol[i],anno$Symbol)
#   if (is.integer(a)&length(a)!=0){
#   res[i] <- a
#   }}
# 
# res1 <- vector(length = nrow(idhuman) )
# 
# for (i in 1: nrow(idhuman)){
#   a<-0
#   print(i)
#   a<- grep(idhuman[i,5],anno$Symbol)
#   if (is.integer(a)&length(a)!=0){
#     res1[i] <- a
#   }}
# x <- anno$HomoloGene.ID[res]
# y <- anno$HomoloGene.ID[res1]
# y[unique(match(x,y))[-1]]
