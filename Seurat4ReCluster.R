rdsFile = "InputFile.rds"       #输入SeuratCluster的rds
outpath = "OutputDir"           #输出路径
#如果输入celllist，一列：取两者交集细胞，两列：取交集并添加初始ident数据
#如果输入selectedClusters，筛选要细分的Cluster
celllist = NULL                 #筛选输入文件中的细胞列表
selectedClusters = ""           #选择要细分的Cluster名称

resolution = 0.8                #分辨率，采用graphbased方法 分群的时候分辨率越大 类型越多 反之越少
minGene = NULL                  #每个细胞最少有n个基因表达，否则过滤掉这个细胞。
#计算MarkerGene参数
MarkerGeneMethod = "wilcox"     #计算markergene的算法
threshold = 0.25                #差异倍数 0.25以下不算显著Marker，会被去除，默认0.25
minpct = 0.1                    #计算MarkerGene表达占总体细胞的最低百分比
onlypos = TRUE                  #只返回positive基因
base=2                          #计算markergene的差异倍数使用2为底

#重降维参数
pcadim = 1:10                   #降维使用的pca维度
computePCs = 50                 #主成分分析返回个数
nneighbors = 30                 #UMAP分析neighbors数量
varGeneNum = 2000               #计算variablegene的数量

library(plyr)
library(reshape2)
library(Seurat)
library(mclust)
library(dplyr)
library(stringr)

setwd(outpath)
seuset <- readRDS(rdsFile);

if(!is.null(celllist)){
	celllistUse = unique(read.table(celllist,sep="\t",header=F)[,1],comment.char = "")
	celllistUse = celllistUse[is.finite(match(celllistUse,colnames(seuset)))]
	seuset = subset(seuset,cells=celllistUse)
}else if(!is.null(selectedClusters)){
    clusters = unlist(strsplit(selectedClusters, "[,]"))
    celllistUse = names(seuset@active.ident[is.finite(match(seuset@active.ident,clusters))])
    seuset = subset(seuset,cells=celllistUse)
}

#对Seurat对象进行重新降维，重跑降维tsne/umap

if(!is.null(minGene)){
    seuset <- subset(x = seuset, subset = nFeature_RNA > minGene)
}
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = varGeneNum)
write.table(VariableFeatures(seuset),file="varGene.txt",sep="\t",row.names=F,quote=F)

scalegene = VariableFeatures(seuset)
model = seuset@commands$ScaleData.RNA$model.use
regressvar = seuset@commands$ScaleData.RNA$vars.to.regress   	

seuset <- ScaleData(object = seuset, features = scalegene, vars.to.regress = regressvar, model.use = model)
seuset <- RunPCA(object = seuset, features = VariableFeatures(seuset), npcs = computePCs, do.print = TRUE, ndims.print = 1:computePCs, nfeatures.print = 10,seed.use=2021)
seuset <- JackStraw(seuset, num.replicate = 100,dims = computePCs)
seuset <- ScoreJackStraw(seuset, dims = 1:computePCs)

seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors, reduction.name = "umap3d", n.components = 3, reduction.key = "umap3d_")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors)
seuset <- RunTSNE(object = seuset,dims = pcadim,reduction.name = "tsne3d",dim.embed=3, reduction.key = "tSNE3d_",check_duplicates = FALSE)
seuset <- RunTSNE(object = seuset,dims = pcadim,check_duplicates = FALSE)

seuset <- FindNeighbors(seuset, reduction = "pca", dims = pcadim)
seuset <- FindClusters(seuset, resolution = resolution)

table = table(Idents(seuset),seuset@meta.data$orig.ident)
table = dcast(data.frame(table),Var1~Var2)
colnames(table)[1]="Cluster"
write.table(table,file=paste0("GraphClust.Statistics.txt"),sep="\t",row.names=F,quot=F)
data = data.frame(Cell = colnames(seuset),Cluster=Idents(seuset))
write.table(data,file=paste0("GraphClust.Summary_Cell.txt"),sep="\t",row.names=F,quot=F)

allmarkers <- FindAllMarkers(object = seuset, only.pos = onlypos, min.pct = minpct, logfc.threshold = threshold, test.use = MarkerGeneMethod,base=base)
write.table(data.frame(allmarkers$gene,allmarkers[,1:6]),file=paste0("GraphClust.AllMarkerGenes.txt"),sep="\t",row.names=F,quot=F)
saveRDS(seuset, file = paste0("GraphClust.SeuratObject.rds"))


