# Laura 1st of March 2021
# script to analise presence/absence of genes from Roary csv file

library(VennDiagram)

# import data and metadata:
dat <- read.csv("./data/gene_presence_absence.csv")
clusters <- read.csv("./data/clusters.csv", as.is = TRUE)

# remove columns that are not very important
dat <- dat[,-5:-14]

# remove hypothetical protein 
dat <- dat[- grep("hypothetical", dat$Annotation ),]
dat <- dat[- grep("hypothetical", dat$Non.unique.Gene.name),]

# remove IS
dat <- dat[- grep("transposase", dat$Annotation),]
dat <- dat[- grep("transposase", dat$Non.unique.Gene.name),]
dat <- dat[- grep("transposase", dat$Gene),]

# produce a matrix for presence/absence of each gene in each strain:
geneMatrix <- dat[ , 5:19]
row.names(geneMatrix) <- dat[ , 1]
geneMatrix <- geneMatrix != ""
# transpose table so that strains are the rows
geneMatrix <- t(geneMatrix)


pdpGeneMatrix <- geneMatrix[clusters[clusters[,2]=="Pdp",1],]
fPddGeneMatrix <- geneMatrix[clusters[clusters[,2]=="fPdd",1],]
tPddGeneMatrix <- geneMatrix[clusters[clusters[,2]=="tPdd",1],]

pangenome_pdp <- colnames(pdpGeneMatrix[,colSums(pdpGeneMatrix)>0])
pangenome_fPdd <- colnames(fPddGeneMatrix[,colSums(fPddGeneMatrix)>0])
pangenome_tPdd <- colnames(tPddGeneMatrix[,colSums(tPddGeneMatrix)>0])

ncol(pdpGeneMatrix)
length(pangenome_pdp)
length(pangenome_fPdd)
length(pangenome_tPdd)

intersect(pangenome_pdp, pangenome_fPdd)

venn.diagram(
  x = list(pangenome_pdp,pangenome_fPdd, pangenome_tPdd),
  category.names = c("Pdp", "fPdd", "tPdd"),
  filename = './plots/venn_diagram.png',
  output = TRUE
)


