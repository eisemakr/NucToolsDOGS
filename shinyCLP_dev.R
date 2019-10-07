
###########################################################################

#This is the CLP shiny App.
#It's pretty cool, try it out

#Just press the green <Run App> Button at the top right side. 

#Pressing <only top 10 species> calculates the variance for every species
#and keeps only the top 10 overall most variant species. 

#In general the kraken data is processed that all rows containing only 0 
#are deleted and then only the top 25 % of all values are used. If 
#<only top 10 species> is set this filtering is not! applied. 

#If I am not saving both plots. Please restart me. 

###########################################################################

if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/00Archive/pacman/pacman_0.4.6.tar.gz",repos=NULL, method="libcurl")
pacman::p_load(lattice, grid, gridExtra, RColorBrewer, viridis, matrixStats,
               reshape2, tidyverse, RMySQL, DT, shiny, factoextra, ggbiplot, limma)
setwd('Y:/IFB_IvD/IB1/iii_Mitarbeiter/Eike Krautter/')

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Import general data

# place the connection here to use the same con for all my queries
con <- dbConnect(MySQL(),user="fungen", password="Test=1234",dbname="CLP", host="10.11.27.118")
rs <- dbSendQuery(con, "SELECT * FROM `Info` WHERE 1;")
#import all CLP sample annotations for selection 
anno <<- fetch(rs)
dbClearResult(rs)
dbDisconnect(con)

anno <- anno[order(anno$sampleID),]
rownames(anno) <- anno$sampleID
anno <- anno[-which(anno$sampleID=='muethe_023'),]
colnames(anno)[c(5,9,10)] <- c('time','depth','library')
  
# load and merge all species_filter_norm files
  
tabl <- 'X:/db02/Sequencing_Projects/Project_chrhar_mouse_ma/final_analysis0.6/countsFiltered/species_filter_norm.txt'
tabl1 <- 'X:/db03/Sequencing_Projects/Project_muethe_mouse_ma/final_analysis0.6.oldDB/countsFiltered/species_filter_norm.txt'
tabl2 <- 'X:/db03/Sequencing_Projects/Project_muethe_mouse_ma_part2/final_analysis0.6/countsFiltered/species_filter_norm.txt'
tabl3 <- 'X:/db02/Sequencing_Projects/Project_walmad_in-vivo_mouse/final_analysis0.6/countsFiltered/species_filter_norm.txt'
  
data <- read.table(tabl,header = T, sep = ",")
data1 <- read.table(tabl1,header = T, sep = ",")
data2 <- read.table(tabl2,header = T, sep = ",")
data3 <- read.table(tabl3,header = T, sep = ",")
  
data <- merge(data,data1,by =0, all = T)
rownames(data) <- data$Row.names
data$Row.names <- NULL
data <- merge(data,data2,by =0, all = T)
rownames(data) <- data$Row.names
data$Row.names <- NULL
data <- merge(data,data3,by =0, all = T)
rownames(data) <- data$Row.names
data$Row.names <- NULL
  
data[is.na(data)] <- 0
data <- data[, order(names(data))]
data <- round(data,2)

### Define the user interface  
ui <- pageWithSidebar(
  # App title ----
  headerPanel("CLP studies"),
  # Sidebar panel for inputs ----
  sidebarPanel(      
    selectInput("variable5", "sampleProject:", multiple = T, c(sort(unique(anno$sampleProject)))),
    selectInput("variable6", "sampleID:", multiple = T, c(unique(anno$sampleID))),
    selectInput("variable1", "mouseID:", multiple = T, c(sort(unique(anno$mouseID)))), 
    selectInput("variable7", "CLPstudy:", multiple = T, c(unique(anno$CLPstudy))), 
    selectInput("variable4", "time:", multiple = T,  c(sort(unique(anno$time)))),
    selectInput("variable2", "sampletype:", multiple = T, c(unique(anno$sampletype))),
    selectInput("variable3", "operation:", multiple = T, c(unique(anno$operation))),
    selectInput("variable8", "antibiotic:", multiple = T, c(unique(anno$antibiotic))),
    sliderInput("variable9", "sequencing depth:", min(anno$depth), max(anno$depth), mean(anno$depth)),
    selectInput("variable10", "library:", multiple = T, c(unique(anno$library))),
    checkboxInput("topten","only top 10 variant species"),
    downloadLink("downloadData", "Download pictures"),
    downloadLink("downloadraw", "Download raw data"),
    width = 2
  ),
  # Main panel for displaying outputs ----
  mainPanel(
    tabsetPanel(
      tabPanel("Selected Samples",dataTableOutput(outputId = 'anno')),
      tabPanel("Species_norm",dataTableOutput(outputId = 'txt')),
      tabPanel("Plot",plotOutput(outputId = 'plot')),
      tabPanel("PCA",plotOutput(outputId = 'pca'))
    )))

### Define the server structure

server <- function(input, output) {
  annoplot = reactive({  
    if (is.null(input$variable1) )  {a <- anno$mouseID}         else{a <- input$variable1}
    if (is.null(input$variable2) )  {b <- anno$sampletype}      else{b <- input$variable2}
    if (is.null(input$variable3) )  {c <- anno$operation}       else{c <- input$variable3}
    if (is.null(input$variable4) )  {d <- anno$time}            else{d <- input$variable4}
    if (is.null(input$variable5) )  {a1 <- anno$sampleProject}  else{a1 <- input$variable5}
    if (is.null(input$variable6) )  {a2 <- anno$sampleID}       else{a2 <- input$variable6}
    if (is.null(input$variable7) )  {a3 <- anno$CLPstudy}       else{a3 <- input$variable7}
    if (is.null(input$variable8) )  {a4 <- anno$antibiotic}     else{a4 <- input$variable8}    
    if (is.null(input$variable9) )  {a5 <- anno$depth}          else{a5 <- input$variable9}
    if (is.null(input$variable10))  {a6 <- anno$library}        else{a6 <- input$variable10}
    
    subset(anno,  (anno$mouseID %in% a ) & (anno$sampletype %in% b ) & (anno$operation %in% c ) & (anno$time %in% d ) & 
                  (anno$sampleProject %in% a1 ) & (anno$sampleID %in% a2 ) & (anno$CLPstudy %in% a3 ) & (anno$antibiotic %in% a4 ) & 
                  (anno$depth <= a5 ) & (anno$library %in% a6 ))
  })
  
  dataplot <- reactive({    
    if (is.null(input$variable1) )  {a <- anno$mouseID}         else{a <- input$variable1}
    if (is.null(input$variable2) )  {b <- anno$sampletype}      else{b <- input$variable2}
    if (is.null(input$variable3) )  {c <- anno$operation}       else{c <- input$variable3}
    if (is.null(input$variable4) )  {d <- anno$time}            else{d <- input$variable4}
    if (is.null(input$variable5) )  {a1 <- anno$sampleProject}  else{a1 <- input$variable5}
    if (is.null(input$variable6) )  {a2 <- anno$sampleID}       else{a2 <- input$variable6}
    if (is.null(input$variable7) )  {a3 <- anno$CLPstudy}       else{a3 <- input$variable7}
    if (is.null(input$variable8) )  {a4 <- anno$antibiotic}     else{a4 <- input$variable8}    
    if (is.null(input$variable9) )  {a5 <- anno$depth}          else{a5 <- input$variable9}
    if (is.null(input$variable10))  {a6 <- anno$library}        else{a6 <- input$variable10}

    dataplottab <- data[,which((anno$mouseID %in% a ) & (anno$sampletype %in% b ) & (anno$operation %in% c ) & 
    (anno$time %in% d )& (anno$sampleProject %in% a1 ) & (anno$sampleID %in% a2 ) & (anno$CLPstudy %in% a3 ) & 
    (anno$antibiotic %in% a4 ) & (anno$depth <= a5 ) & (anno$library %in% a6 ))]
    
    dataplottab <- dataplottab[which(rowSums(dataplottab)>0), ]
    if (input$topten==T){
      dataplottab <- cbind(dataplottab, var = rowVars(as.matrix(dataplottab)))
      dataplottab <- dataplottab[order(dataplottab$var, decreasing = T),]
      dataplottab <- dataplottab[c(1:10),]
      dataplottab$var <- NULL
      dataplottab
    }else{
      dataplottab
    }
    })

  
  
  plotspecies <- reactive({
    
    dataplot1 <- dataplot()
    dataplot1 <- cbind(dataplot1,rownames(dataplot1))
    dataplot1 <- melt(dataplot1)
    colnames(dataplot1) <- c('bacteria','variable','value') 
    dataplot1 <- dataplot1[-which(dataplot1$value == 0),]
    if (input$topten==F){
      dataplot1 <- dataplot1[-which(dataplot1$value < quantile(dataplot1$value)[3]),]
    }else{
    }
    if (length(unique(dataplot1$bacteria ))<= 74){
      ggplot(data=dataplot1)+
        geom_bar(aes(x=variable,y=value,fill=bacteria), stat = 'identity')+
        scale_fill_manual(values=col_vector)+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }else{
      ggplot(data=dataplot1)+
        geom_bar(aes(x=variable,y=value,fill=bacteria), stat = 'identity')+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        theme(legend.text = element_text(colour="blue", size=8, face="bold"))+  
        guides(fill=guide_legend(ncol = 2))  
    }
  })

  plotpca <- reactive({
    
    pca <- prcomp(t(dataplot()),scale=FALSE,retx=TRUE)
    pca_use <- as.data.frame(pca$x)
    pca_use <- cbind(pca_use,as.factor(annoplot()$sampletype),as.factor(annoplot()$time),as.factor(sub(' ','',annoplot()$antibiotic)))
    colnames(pca_use)[(ncol(pca_use)-2):ncol(pca_use)] <- c('source', 'time', 'AB')
    rownames(pca$x) <- sub('chrhar_','c',rownames(pca$x))
    rownames(pca$x) <- sub('muethe_','m',rownames(pca$x))
    rownames(pca$x) <- sub('walmad_','w',rownames(pca$x))
    
    plot1 <- fviz_eig(pca)
    plot2 <- ggplot(pca_use)+
      geom_text(aes(x=PC1,y=PC2, colour = time, size = AB), label=rownames(pca$x))+
      theme_classic()
    
    plot3 <- ggplot(pca_use)+
      geom_text(aes(x=PC1,y=PC3, colour = time, size = AB), label=rownames(pca$x))+
      theme_classic()
    
    plot4 <- ggplot(pca_use)+
      geom_text(aes(x=PC2,y=PC3, colour = time, size = AB), label=rownames(pca$x))+
      theme_classic()
    plot4
    grid.arrange(plot1, plot2, plot3, plot4, ncol=2)
    
  })
  
  plotpdf <- reactive({
    list(p1 =  plotspecies(), p2 = plotpca())
  })
  
  output$pca <- renderPlot({plotpca()},height = 1200,width = 1800)
  output$plot <- renderPlot({plotspecies()},height = 1200,width = 1800)
  output$anno = DT::renderDataTable({annoplot()})
  output$txt <- DT::renderDataTable({dataplot()})
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if (input$topten == T){
        top <- 'top10'
      }else{
        top <- 'allspec'
      }
      gsub("not","",gsub("-","",paste(input$variable1,input$variable2,input$variable3,input$variable4,
            input$variable11,input$variable21,input$variable31,input$variable41,top,".pdf", sep = ""),fixed=T),fixed=T)
    },
    content = function(file) {
      pdf(file, width = 12)
      invisible(lapply(plotpdf(), print))
      dev.off()
    }
  )
  
  output$downloadraw <- downloadHandler(
    filename = function() {
      if (input$topten == T){
        top <- 'top10'
      }else{
        top <- 'allspec'
      }
      gsub("not","",gsub("-","",paste(input$variable1,input$variable2,input$variable3,input$variable4,
                                      input$variable11,input$variable21,input$variable31,input$variable41,top,"_rawdata.txt", sep = ""),fixed=T),fixed=T)
    },
    content = function(file) {
      write.table(dataplot(),file)

    }
  )
}

### start the App 

app <- shinyApp(ui,server)
runApp(app, launch.browser = TRUE)

