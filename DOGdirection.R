#This script takes the output of pos2reg (.bed file with DOGS) and the output of DOGS2Genes
#(A big bed file with intersection results) the direction of the DOG stored in the first output is 
#added to the second file and a new output is produced. 

if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.4.6.tar.gz",repos=NULL, method="libcurl")
pacman::p_load( biomaRt, dplyr )

args <- commandArgs(TRUE)
p2data <- as.character(args[1])   # path to DOG file
orga <- as.character(args[2])   # path to Gene file

if(is.na(args[1])|is.na(args[2])){
  print('please specify both arguments!')
  quit()
}
