library(biomaRt)
anno <- read.table('X:/db05/Eike/mouse_human-txt.txt', fill = T, sep = '\t')
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
