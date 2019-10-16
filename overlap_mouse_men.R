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

#for development

p2datahuman <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/Control/con1_con2/DOGS2dirgene1000.bed'
p2datamouse <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/cutoff30/DOGS2dirgene1000.bed'
wd <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/Control/dev'



if(is.na(args[1])|is.na(args[2])|is.na(args[3])){
  print('please specify all 3 arguments!')
  quit()
}

setwd(wd)

if (getwd()!= wd) {
  print('working directory could not be set. Check pls thats essential!')
  quit()
}

datahuman <- read.table(p2datahuman)
datamouse <- read.table(p2datamouse)

print('using grch37!!! Change my code if you need the latest human genome.')
ensembl.human = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        host="grch37.ensembl.org",
                        path="/biomart/martservice",
                        dataset="hsapiens_gene_ensembl")
idhuman <- getBM(attributes = c('mmusculus_homolog_associated_gene_name', 'mmusculus_homolog_chromosome'), 
                 filters = 'ensembl_transcript_id',
                 values = datahuman$V4,
                 mart = ensembl.human)
ensembl.mouse= useMart("ensembl",dataset="mmusculus_gene_ensembl")
idmouse <- getBM(attributes = c('mgi_symbol','chromosome_name','ensembl_transcript_id'), 
                 filters = 'ensembl_transcript_id',
                 values = datamouse$V4,
                 mart = ensembl.mouse)



anno <- read.table('X:/db05/Eike/mouse_human-txt.txt', fill = T, sep = '\t')

ensembl <- useMart('ensembl')
ensembl.human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
ensembl.mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
      filters = 'hgnc_symbol',
      values = c('TP53', 'XPC', 'DDB1', 'APEX1'),
      mart = ensembl.human)
getBM(attributes = c('mmusculus_homolog_associated_gene_name', 'mmusculus_homolog_chromosome'),
      filters = 'hgnc_symbol',
      values = c('E2F3'),
      mart = ensembl.human)
my.mouse.genes <- c('Frem2', 'Kmt2d', 'Scn8a', 'Abcg1', 'Acvr1b', 'Flnc', 'Lama2', 'Myh7b', 'Myo10', 'Ryr3')
getLDS(attributes = c('mgi_symbol', 'chromosome_name'),
       filters = 'mgi_symbol', values = my.mouse.genes, mart = ensembl.mouse,
       attributesL = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position') , martL = ensembl.human)


quit()

















data <- read.table('X:/db05/Eike/genesmuethe_mirsi.txt', fill = T, sep = '\t')

res <- vector(length = length(data$V3) )

for (i in 1: length(data$V3)){
  a<-0
  print(i)
  a<- grep(data$V3[i],anno$V11)
  if (is.integer(a)&length(a)!=0){
  res[i] <- a
  }}

res1 <- vector(length = length(data$V2) )

for (i in 1: length(data$V2)){
  a<-0
  print(i)
  a<- grep(data$V2[i],anno$V11)
  if (is.integer(a)&length(a)!=0){
    res1[i] <- a
  }}
x <- anno$V1[res]
y <- anno$V1[res1]
y[unique(match(x,y))]
