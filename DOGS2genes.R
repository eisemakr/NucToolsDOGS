#This script takes a list of DOGS in .bed-file format and intersects it with the corresponding transcription start sites extended to a size defined
#by size parameter. The list of TSS is automatically downloaded by biomaRt from ENSEMBL. Further stats are printed and intersections with regulatory 
#features from ENSEMBL are down but only their stats were printed. All IDs are ENSEMBL transcript IDs. This script was mainly made due to the big 
#differences in IDs from different manually downloaded files and the imposibillity to convert those. 
#Working directory is always target directory. Please specify full path to target file. output file will be in directory of target file. 

#################################################################################

#Load packages, check inputs and print inputs

#################################################################################

if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.4.6.tar.gz",repos=NULL, method="libcurl")
pacman::p_load( biomaRt, dplyr )
# library(dplyr)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biobase")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
# 
# library(biomaRt)


args <- commandArgs(TRUE)
p2data <- as.character(args[1])   # path to target file
orga <- as.character(args[2])   # eighter mouse or human 
size <- as.numeric(args[3])   # size to which transcription start sites will be extended

if(is.na(args[1])|is.na(args[2])|is.na(args[3])){
  print('please specify all 3 arguments!')
  quit()
}

#################################################################################

#Check inputs, create respective mart object, remove all unusfull chromosomes, 
#extend by size, save retrieved data in .bed format
#not implemented -- check if bedtools is installed

#################################################################################

if (size>2000|size<300){   
  print('size is out of recommended range. Please reconsider your choice!')
  quit()
}

wd <- paste(sub('\\/\\w*\\.\\w*', '', p2data, perl=T), sep = '')
print('Working directory is:')
setwd(wd)
print(getwd())

if (getwd()!= wd) {
  print('working directory could not be set. A regex is used to extract it eigther your path is wrong or my regex bad. Check both pls.')
  quit()
}

if (orga=='mouse'){
   print('get genelist')
   ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
   transreg <- getBM(attributes = c('chromosome_name','strand', 'transcription_start_site','ensembl_transcript_id','transcript_length'), mart = ensembl)
   transreg <- transreg[which(transreg$chromosome_name %in% c(1:19, 'X', 'Y')),]
   print('get exonlist') # I never tested the exon part for mouse
   exons <- getBM(attributes = c('chromosome_name','exon_chrom_start', 'exon_chrom_end', 'ensembl_transcript_id'), mart = ensembl)
   exons <- exons[which(exons$chromosome_name %in% c(1:19, 'X', 'Y')),]
   print('get regfeatures')
   ensi <- useMart("ENSEMBL_MART_FUNCGEN", dataset = 'mmusculus_regulatory_feature')
   reganno <- getBM(attributes = c('chromosome_name','chromosome_start', 'chromosome_end',"feature_type_name"), mart = ensi)
   #reganno <- reganno[which(reganno$chromosome_name %in% c(1:19, 'X', 'Y')),]# I never tested this but actually it would be useful to cut away everything we don't need

 }else if (orga=='human'){
   print('get genelist')
   #This is for hg19!!!!! The code has to be changed if we use GRCh38
   ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
   transreg <- getBM(attributes = c('chromosome_name','strand', 'transcription_start_site','ensembl_transcript_id','transcript_length'), mart = ensembl)
   transreg <- transreg[which(transreg$chromosome_name %in% c(1:22, 'X', 'Y')),]
   #dele <- getBM(attributes = c('chromosome_name','transcript_start', 'transcript_end','ensembl_transcript_id'), mart = ensembl)
   #dele <- dele[which(dele$chromosome_name %in% c(1:22, 'X', 'Y')),]
   print('get exonlist')
   exons <- getBM(attributes = c('chromosome_name','exon_chrom_start', 'exon_chrom_end', 'ensembl_transcript_id'), mart = ensembl)
   exons <- exons[which(exons$chromosome_name %in% c(1:22, 'X', 'Y')),]
   print('get regfeatures')
   ensi <- useMart("ENSEMBL_MART_FUNCGEN", host="grch37.ensembl.org", dataset = 'hsapiens_regulatory_feature')
   reganno <- getBM(attributes = c('chromosome_name','chromosome_start', 'chromosome_end', 'feature_type_name'), mart = ensi)
   #reganno <- reganno[which(reganno$chromosome_name %in% c(1:22, 'X', 'Y')),]# I never tested this but actually it would be useful to cut away everything we don't need

 }else{
   print('this message should never occur!!! you specified an incorrect orga but it wasnt filtered away')
   quit()
 }

#Using the transcription starts, strandedness and size parameter to get the regions upstream the transcription start in .bed format

transreg1 <- subset(transreg,transreg$strand==1)
transcripts1 <- as.data.frame(cbind(transreg1$chromosome_name,transreg1$transcription_start_site,
                    transreg1$transcription_start_site + transreg1$transcript_length, transreg1$ensembl_transcript_id),
                    stringsAsFactors = F)
transreg1$strand <- transreg1$transcription_start_site - size
transreg2 <- subset(transreg,transreg$strand==-1)
transcripts2 <- as.data.frame(cbind(transreg2$chromosome_name, transreg2$transcription_start_site - transreg2$transcript_length,
                                 transreg2$transcription_start_site, transreg2$ensembl_transcript_id),
                             stringsAsFactors = F)
transreg2$strand <- transreg2$transcription_start_site
transreg2$transcription_start_site <- transreg2$transcription_start_site + size

transreg <- rbind(transreg1,transreg2)
transreg <- transreg[order(transreg$chromosome_name,transreg$strand),]

transcripts <- rbind(transcripts1,transcripts2)
transcripts <- transcripts[order(transcripts$V1,transcripts$V2),]

#################################################################################

#Read the file with the DOGS, convert it to small .bed and write all annotation
#files to the directory in order to make the report reproducible.

#################################################################################

data <- read.table(p2data, fill = T)
data <- data[,1:4]
data[,1] <- sub('chr', '', data[,1])

annoname <- paste(Sys.Date(), size, 'upstreamGenAnnotation.bed', sep = '')
transcriptsname <- paste(Sys.Date(), 'GenAnnotation.bed', sep = '')
regname <- paste(Sys.Date(), 'RegAnnotation.bed', sep = '')
dataname <- 'DOGSsmallbed.bed'
exonname <- 'exonshg19.bed'

write.table(transreg , file = annoname, sep = '\t', quote = F, row.names = F, col.names = F)
write.table(reganno , file = regname, sep = '\t', quote = F, row.names = F, col.names = F)
write.table(data , file = dataname, sep = '\t', quote = F, row.names = F, col.names = F)
write.table(exons , file = exonname, sep = '\t', quote = F, row.names = F, col.names = F)
write.table(transcripts , file = transcriptsname, sep = '\t', quote = F, row.names = F, col.names = F)

#################################################################################

#run bedtools intersect from the shell to create a file of DOGS and genes

#################################################################################

print(paste('bedtools intersect -wo -a ',annoname, ' -b ',dataname,' > DOGS2dirgene',size,'.bed', sep = ''))
print(paste('bedtools intersect -wo -a ',regname, ' -b ',dataname,' > DOGS2regs.bed', sep = ''))

system(paste('bedtools intersect -wo -a ',annoname, ' -b ',dataname,' > DOGS2dirgene',size,'.bed', sep = ''))
system(paste('bedtools intersect -wo -a ',regname, ' -b ',dataname,' > DOGS2regs.bed', sep = ''))

nbDOGS2Genes <- system(paste('bedtools intersect -wo -a ',annoname, ' -b ',dataname,' | wc -l ', sep = ''),
                      intern = T)
nbDOGS2Transcripts <- system(paste('bedtools intersect -wo -a ',transcriptsname , ' -b ',dataname,' | wc -l ', sep = ''),
                       intern = T)
nbDOGS2exon <- system(paste('bedtools intersect -b ', dataname,' -a ', exonname,' | wc -l ', sep = ''), intern = T)

#receive information from the intersection with the regulatory features file 

nbDOGS2Pro <- system('cat DOGS2regs.bed | grep Promoter | grep -v "Promoter Flanking Region" | wc -l ', intern = T)
nbDOGS2Profl <- system('cat DOGS2regs.bed | grep "Promoter Flanking Region" | wc -l ', intern = T)
nbDOGS2Enh <- system('cat DOGS2regs.bed | grep "Enhancer" | wc -l ', intern = T)
nbDOGS2CTCF <- system('cat DOGS2regs.bed | grep "CTCF Binding Site" | wc -l ', intern = T)
nbDOGS2TF <- system('cat DOGS2regs.bed | grep "TF binding site" | wc -l ', intern = T)
nbDOGS2OC <- system('cat DOGS2regs.bed | grep "Open chromatin" | wc -l ', intern = T)

#Put the results to a textfile. The textfile will not be overwritten by running the program again. 

txt <- c(paste('number of DOGS ', size ,' bases upstream of genes is: ', nbDOGS2Genes, sep = ''),
  paste('number of DOGS inside annotated promoters is: ', nbDOGS2Pro, sep = ''),
  paste('number of DOGS inside annotated promoter flanking regions is: ', nbDOGS2Profl, sep = ''),
  paste('number of DOGS inside enhancers is: ', nbDOGS2Enh, sep = ''),
  paste('number of DOGS inside CTCF binding sites is: ', nbDOGS2CTCF, sep = ''),
  paste('number of DOGS inside TF binding sites is: ', nbDOGS2TF, sep = ''),
  paste('number of DOGS inside open chromatin is: ', nbDOGS2OC, sep = ''),
  paste('number of DOGS inside transcripts is: ', nbDOGS2Transcripts, sep = ''),
  paste('number of DOGS inside exons is: ', nbDOGS2exon, sep = ''),
  paste(''),
  paste('All numbers were derived by using standard settings of bedtools intersect (1 base overlap?)'))

print(txt)

sink(file = 'intersectionstats.txt', append = T)
print(Sys.Date())
print(Sys.time())
print('inputs were:')
print(paste(as.character(args[1]), as.character(args[2]), as.numeric(args[3])))
print(txt)
sink()

quit()
