#This script takes the output of pos2reg (.bed file with DOGS) and the output of DOGS2Genes
#(A big bed file with intersection results) the direction of the DOG stored in the first output is 
#added to the second file and a new output is produced. 

if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.4.6.tar.gz",repos=NULL, method="libcurl")
pacman::p_load( dplyr )

args <- commandArgs(TRUE)
p2DOG <- as.character(args[1])   # path to DOG file
p2gene <- as.character(args[2])   # path to Gene file

#################################################################################
# For development
p2DOG <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/Control/dev/con1_con2/TestDOGS.bed' 
p2gene <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/Control/dev/con1_con2/DOGS2dirgene1000.bed'
#################################################################################

if(is.na(args[1])|is.na(args[2])){
  print('please specify both arguments!')
  quit()
}

DOGS <- read.table(p2DOG, fill = T, header = F)
GENES <- read.table(p2gene, fill = T, header = F)

#################################################################################
#Extract the DOG identifiers from the gene file
#################################################################################

DOGSlist <- as.character(unique(GENES$V9))

#################################################################################
#Subset the DOG file by only the DOGS matching a gene 
#################################################################################

DOGSwGENES <- DOGS[which(DOGS$V4 %in% DOGSlist),]

DOGDIRobject <- merge(GENES,DOGSwGENES, by.x = 'V9', by.y = 'V4', sort = F)

outnm <- sub('.bed','',p2gene)
outnm <- paste(outnm, '_dir.bed', sep = '')

write.table(DOGDIRobject, outnm,  quote = F, row.names = F, col.names = F)

quit()
