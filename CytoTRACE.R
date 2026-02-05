library(Seurat)
library(CytoTRACE)
library(ggplot2)
library(cowplot)

rdsFile = "cluster.rds"

seuset=readRDS(rdsFile);
data = GetAssayData(seuset, slot = "counts")
data = as.matrix(data[,colnames(seuset)])
results <- CytoTRACE(data,ncores=1)
CytoValue = results$CytoTRACE
Cytodata = data.frame(Cell = names(CytoValue),CytoTRACE=CytoValue)
seuset = AddMetaData(seuset,Cytodata)
CytoGene = results$cytoGenes
CytoGene = sort(CytoGene,decreasing = T)
Genedata = data.frame(Gene = names(CytoGene),GenesCorrelated=CytoGene)

forplot = data.frame(seuset@meta.data$CytoTRACE,seuset@active.ident)
colnames(forplot) = c("CytoTRACE","Cluster")
gg = ggplot(data = forplot, mapping = aes(x= Cluster,y = CytoTRACE, col=Cluster,fill=Cluster))+geom_boxplot(alpha=0.5,outlier.size = -1)

plotgeneNum = 10
forplot = rbind(head(Genedata,plotgeneNum),tail(Genedata,plotgeneNum))
forplot$TF = forplot$GenesCorrelated<=0
gg = ggplot(data=forplot,mapping = aes(y=GenesCorrelated,x=Gene,fill=TF))+geom_col()+coord_flip()
