# if (!requireNamespace("BiocManager", quietly = TRUE))
#   
# install.packages("BiocManager")
# BiocManager::install("KEGGprofile")

library(biomaRt)
library(KEGGREST)
library(KEGGprofile)
#download_KEGGfile(pathway_id="all",species='hsa')

data(pro_pho_expr)
data(pho_sites_count)
ls()
pro_pho_expr[1:3,1:4]

data <- read.table('W:/db05/Eike/mirsi/Analysis/stabfuz/Control/con1_con2/DOGS2dirgene1000.bed')
data$V4
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
nms <- getBM(attributes = c('ensembl_transcript_id','entrezgene_id'),
                  filters = 'ensembl_transcript_id',
                  values = data$V4,
                  mart = ensembl)

genes <- unique(nms$entrezgene_id)
genes<- as.character(genes[-1])
res<-find_enriched_pathway(genes,species='hsa')
res[[1]][,c(1,5)]


data <- read.table('W:/db05/Eike/mirsi/Analysis/stabfuz/First/con0_cas0/DOGS2dirgene1000.bed')
data$V4
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
nms <- getBM(attributes = c('ensembl_transcript_id','entrezgene_id'),
             filters = 'ensembl_transcript_id',
             values = data$V4,
             mart = ensembl)

genes <- unique(nms$entrezgene_id)
genes<- as.character(genes[-1])
res<-find_enriched_pathway(genes,species='hsa')
res[[1]][,c(1,5)]

data <- read.table('W:/db05/Eike/mirsi/Analysis/stabfuz/cutoff30/DOGS2dirgene1000.bed')
data$V4
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
nms <- getBM(attributes = c('ensembl_transcript_id','entrezgene_id'),
             filters = 'ensembl_transcript_id',
             values = data$V4,
             mart = ensembl)
genes <- unique(nms$entrezgene_id)
genes<- as.character(genes[-1])
res<-find_enriched_pathway(genes,species='mmu')
res[[1]][,c(1,5)]
geneList <- getBM(attributes = c('entrezgene_id'),
                  #filters = 'ensembl_transcript_id',
                  #values = data$V4,
                  mart = ensembl)
all <- factor(as.integer (geneList$entrezgene_id %in% genes))
names(all) <- as.vector(geneList$entrezgene_id)

# require(DOSE)
# require(clusterProfiler)
# data(geneList)
# gene = names(geneList)[abs(geneList) > 2]
#          enrichKEGG(gene = genes,
#                     organism = "mmu",
#                     keyType = "kegg",
#                     pvalueCutoff = 0.01,
#                     pAdjustMethod = "BH",
#                     qvalueCutoff = 0.05,
#                     minGSSize = 5)

#          
# GO <- enrichGO(gene = genes, 
#         ont = "BP",
#         universe = as.character(geneList$entrezgene_id),
#         OrgDb = org.Mm.eg.db,
#         pAdjustMethod = "BH",
#         pvalueCutoff  = 0.01,
#         qvalueCutoff  = 0.05,
#         minGSSize = 2
#         )
# 
# summary(GO)
# 
# getFnAnot_genome(geneList = genes,
#                  email = 'eike.krautter@igb.fraunhofer.de',
#                  idType = "ENTREZ_GENE_ID",
#                  listName = "auto_list",
#                  count = 1L,
#                  PVal = 1,
#                  getKEGG = FALSE)
