input = "Input_Dir"      #输入文件，CellRanger的矩阵路径
Doubletrate = 0.1     	 #双细胞率,手动调整
outpath = "Output_Dir"   #输出路径

library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)
library(getopt)
library(stringr)
dir.create(outpath)
setwd(outpath)
#载入表达矩阵
allcounts = Read10X(data.dir=input)
#计算双细胞
seuset <- CreateSeuratObject(counts = allcounts,min.cells = 3, min.features = 0)
seuset <- NormalizeData(seuset,normalization.method = "LogNormalize", scale.factor = 10000)
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = 2000)
scalegene = VariableFeatures(seuset)
seuset <- ScaleData(object = seuset,features = scalegene)
sct = F
pc.num = 1:10
nExpadj="true"
resolution = 0.8
pN = 0.25
seuset <- RunPCA(seuset, verbose = F)
seuset <- FindNeighbors(seuset, dims = pc.num) %>% FindClusters(resolution = resolution)
sweep.res.list <- paramSweep_v3(seuset, PCs = pc.num, sct = sct)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
homotypic.prop <- modelHomotypic(seuset$seurat_clusters)
nExp_poi <- round(Doubletrate*ncol(seuset))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp = nExp_poi.adj

seuset <- doubletFinder_v3(seuset, PCs = pc.num, pN = pN, pK = pK_bcmvn, nExp = nExp, reuse.pANN = F, sct = sct)
valuename = paste0("pANN_",pN,'_',pK_bcmvn,'_',nExp)
groupname = paste0("DF.classifications_",pN,'_',pK_bcmvn,'_',nExp)
doubletres = cbind(rownames(seuset@meta.data),seuset@meta.data[,valuename],seuset@meta.data[,groupname])
colnames(doubletres) <- c("Cell","DoubletValue@value","Group@cluster")
write.table(doubletres,paste0("DoubletFinder.result.txt"),sep="\t",quote=F,row.names=F)
