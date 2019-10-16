# This script is a complete follow-up of the comp2cond script by Evgeniy. 
# It takes the output, applies a refiltering for regions which are unoccupied in one condition
# converts positions into regions and calculates the average <<relative error>> of each region. 
# We would call regions DOGS. The outputis a gzip .bed-file and will be printed to the same directory. 
# Take care that the output of comp2cond.pl is sometimes missformatted due to missing chromosome ID. 
# This script is build to be run with the -wE flag from comp2cond. For every DOG the origin file name
# is reported so its "direction" can be traced back. 

#################################################################################

#Load packages, check inputs and print inputs

#################################################################################

if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.4.6.tar.gz",repos=NULL, method="libcurl")
pacman::p_load( stringr )

args <- commandArgs(TRUE)
wd <- args[1]   # working directory, with only gziped files from comp2cond. The scirpt greps for txt.gz.
t <- as.numeric(args[2])    # threshold for region size. Smaller regions will be discarded. (I recommend 30)
g <- as.numeric(args[3])    # maximal gap size for regions to be merged
f <- as.numeric(args[4])    # filtering threshold of minimal mean value of average profile. Where one average profile is 0.
                # This has to be tested somehow, because it can be different due to different nature of data. 
out <- args[5]  # A name for the final output file.

if(is.na(args[1])|is.na(args[2])|is.na(args[3])|is.na(args[4])|is.na(args[5])){
  print('please specify all 5 arguments!')
  quit()
}

#################################################################################
#For testing use these variables
#################################################################################

# t <- 30
# g <- 10
# f <- 0.2
# out <- 'testdelte.bed'
# wd <- 'W:/db05/Eike/mirsi/Analysis/stabfuz/First/con0_cas0'
# i <- 21

#################################################################################

print('arguments are:')
print(paste('working directory is: ',wd))
print('Please make sure that the working directory contains only gzip files to be evaluated.')
print(paste('minimal region size is: ',t))
print(' 30 is recommended.')
print(paste('maximal gap size between regions is: ',g))
print(' I mostly used 10.')
print(paste('minimal mean of average profile in unoccupied sites is: ',f))
print(' 0.5 is a good start. But this highly depends on the data used.')
print(paste('name for the .bed-file containing all regions is: ',out))
print(' There are no recommendations for naming your files :) But it should end with .bed')

setwd(wd)

#################################################################################

#Define the function to convert position-wise format to region-wise format. 

#################################################################################

pos2reg <- function(x,t,g){

if (missing(g)) {
    message("no gap threshold defined. The threshold of merging regions is set to 10 automatically!")
    g <- 10
}  
  
x <- x[,3]
d <- diff(x)
d[length(d)+1] <- NA

x<-as.data.frame(cbind(x,d))

endpos <- x[which(x$d!=1),]$x+1 # end positions

size <- diff(as.numeric(rownames(x[which(x$d!=1),]))) 
size <- c(as.numeric(rownames(x[which(x$d!=1),]))[1],size) # sizevector

positions <- cbind(endpos-(size+1),endpos,size,c(NA,endpos[-length(endpos)]),c(NA,(endpos-size)[-length(endpos)]))
positions <- cbind(positions,positions[,1]-positions[,4])
positions <- positions[-c(which(positions[,6]<g)-1),] 
positions[,1][which(positions[,6]<g)] <- positions[,5][which(positions[,6]<g)]
positions[,3][which(positions[,6]<g)] <- positions[,2][which(positions[,6]<g)]-positions[,1][which(positions[,6]<g)]
positions <- positions[,c(1:3)]
colnames(positions) <- c("start","end","size")

if (missing(t)) {
  message("no threshold defined, all regions will be printed!")
  message(paste("number of regions is", nrow(positions)))
  return(positions)
  stop()
}
message(paste("number of regions removed due to low size is:",length(which(positions[,3]<t))))
positions <- subset(positions, positions[,3]>=t)
message(paste("all regions smaller or equal to",t,"will be removed!"))
message(paste("number of regions is:", nrow(positions)))
return(positions)

}

#################################################################################

#Acctual execution of all steps, looped for all files in a folder

#################################################################################

files <- list.files(pattern = "txt.gz")
nms <- sub('.txt.gz','',files)
#nms <- sub('chr\\d\\d|chr\\d|chrX','',nms)
print('detected files for follow-up analysis are:')
print(files)

for (i in 1:length(files)) {
  print(i)
  data <- read.table(files[i])

  ## filtering

  data1 <- data[which((data$V4==1 & data$V5>f) | (data$V4<1 & data$V4>0)),]
  data2 <- data[which((data$V4==-1 & data$V8>f) | (data$V4>-1 & data$V4<0)),]
  datanew <- rbind(data1, data2)

  ## Applying pos2reg

  reg <- pos2reg(datanew,t,g)
  reg <- as.data.frame(reg)
  
  ## Extracting the average <<rel.error>> for every region.
  # The average relerr of a DOG is the relerr of the experimental being higher in mean occupancy. 
  # If there would be 2 DOGS where high/low would alterate consecutively, we would get a mess.
  
  dir1 <- which(data$V5>data$V8)
  dir2 <- which(data$V5<data$V8)
  
  datacond <- rbind(data[dir1,1:7],data[dir2,c(1:4,8:10)])
  colnames(datacond) <- c('V1','V2','V3','V4','V5','V6','V7')
  
  if (nrow(reg)!=0){
    for (a in 1:nrow(reg)){
      reg$avrelerr[a] <- mean(subset(datacond,datacond$V2>reg$start[a] & datacond$V2<reg$end[a])$V7)
      reg$sdrelerr[a] <- sd(subset(datacond,datacond$V2>reg$start[a] & datacond$V2<reg$end[a])$V7)
    }
    reg <- cbind(chr=str_match(files[i],"chr\\d\\d|chr\\d|chrX"),reg,files[i])
    write.table(reg,file=gzfile(paste(sub(".txt.gz","",files[i]),"_reg",".bed.gz", sep="")), quote = F, col.names = F, sep = '\t', row.names = F )
  }else{
    print('Nothing to be merged. Rare case but happens from time to time.')}
}

#################################################################################

#Creation of the final .bed file

#################################################################################

system('zcat *.bed.gz > temp.bed')
totalbed <- read.table('temp.bed')
totalbed[,8] <- totalbed[,4] 
totalbed[,4] <- paste('DOG', 1:nrow(totalbed), sep="")
totalbed <- totalbed[,c(1,2,3,4,8,5,6,7)]

write.table(totalbed, file = out, quote = F, col.names = F, sep = '\t', row.names = F)
system('rm temp.bed')
quit()
