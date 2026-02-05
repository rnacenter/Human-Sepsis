rdsFile = "InputFile.rds"       #输入SeuratCluster的rds
outpath = "OutputDir"           #输出路径
#如果输入celllist，一列：取两者交集细胞，两列：取交集并添加初始ident数据
#如果输入selectedClusters，筛选要细分的Cluster
celllist = NULL                 #筛选输入文件中的细胞列表
selectedClusters = ""           #选择要细分的Cluster名称
resolution = 0.8                #分辨率，采用graphbased方法 分群的时候分辨率越大 类型越多 反之越少
minGene = NULL                  #每个细胞最少有n个基因表达，否则过滤掉这个细胞。
minCell = 0						#每个样本至少有n个细胞才纳入计算
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


library(Seurat)
library(reshape2)
library(dplyr)
library(stringr)
library(mclust)
library(SingleCellExperiment)
library(scater)
library(plyr)
library(scran)
library(cowplot)
library(batchelor)

dir.create(outpath)
setwd(outpath)

seuset= readRDS(rdsFile)
regressvar = seuset@commands$ScaleData.RNA$vars.to.regress

if(!is.null(minGene)){
  seuset <- subset(x = seuset, subset = nFeature_RNA > minGene)
}

if(!is.null(celllist)){
  #以用户提供的细胞列表，去获取seuratrds中的交集细胞
	celllist = read.table(celllist,header=T,sep="\t")
	rownames(celllist) = celllist[,1]
	celllistUse = unique(celllist[,1])
	celllistUse = celllistUse[is.finite(match(celllistUse,colnames(seuset)))]
	seuset = subset(seuset,cells=celllistUse)
  	if(ncol(celllist)>1){
    	realcelllist = celllist[celllistUse,]
  	}
}else if(!is.null(selectedClusters)){
	#以用户提供的cluster，去获取seuratrds中的交集细胞
	clusters = unlist(strsplit(selectedClusters, "[,]"))
	celllistUse = names(seuset@active.ident[is.finite(match(seuset@active.ident,clusters))])
	seuset = subset(seuset,cells=celllistUse)
}


#将seuset按照样本拆分成list
seuset.list <- SplitObject(seuset, split.by = "orig.ident")

cellNum = lapply(seuset.list, ncol)
#保留符合最小细胞数要求的seurat对象
seuset.list = seuset.list[cellNum>minCell]

#合并样本
seuset = merge(
    x = seuset.list[[1]],
    y = seuset.list[2:length(x = seuset.list)]
)
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = varGeneNum)
scalegene = VariableFeatures(seuset)
seuset <- ScaleData(object = seuset, features = scalegene, vars.to.regress = regressvar, model.use = model)
write.table(VariableFeatures(seuset),file="varGene.txt",sep="\t",row.names=F,quote=F)

#把seurat对象list转成SingleCellExperiment对象的list
objects.sce <- lapply(
    X = seuset.list,
    FUN = function(x, f) { 
      return(as.SingleCellExperiment(x = subset(x = x, features = VariableFeatures(seuset))))
    },
    f = VariableFeatures(seuset)
)

#开始进行MNN标化
print("RunMNN")
out <- do.call(fastMNN, c(objects.sce, list(k = kValue, d = computePCs)));

outcorrected = reducedDim(out)
rownames(outcorrected)=colnames(seuset)
colnames(outcorrected)=paste0("MNN_",1:ncol(outcorrected))
featureloading = as.matrix(rowData(out))
colnames(featureloading)=paste0("MNN_",1:ncol(featureloading))
#在seuset添加mnn结果
seuset[["mnn"]] <- CreateDimReducObject(
    embeddings = outcorrected,
    loadings = featureloading,
    assay = DefaultAssay(object = seuset),
    key = "mnn_"
)
write.table(data.frame(CellName=rownames(seuset@reductions$mnn@cell.embeddings),seuset@reductions$mnn@cell.embeddings),file="mnn.txt",sep="\t",quote=F,row.names=F)

metadata = c("orig.ident")

print("RunUMAP")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors, reduction.name = "umap3d", n.components = 3, reduction.key = "umap3d_",reduction="mnn")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors,reduction="mnn")


write.table(data.frame(CellName=rownames(seuset@reductions$umap@cell.embeddings),seuset@reductions$umap@cell.embeddings),file="umap.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$umap3d@cell.embeddings),seuset@reductions$umap3d@cell.embeddings),file="umap3d.txt",sep="\t",quote=F,row.names=F)


print("RunTSNE")
seuset <- RunTSNE(object = seuset,dims = pcadim,reduction.name = "tsne3d",dim.embed=3, reduction.key = "tSNE3d_",check_duplicates = FALSE,reduction="mnn")
seuset <- RunTSNE(object = seuset,dims = pcadim,check_duplicates = FALSE,reduction="mnn")

write.table(data.frame(CellName=rownames(seuset@reductions$tsne@cell.embeddings),seuset@reductions$tsne@cell.embeddings),file="tsne.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$tsne3d@cell.embeddings),seuset@reductions$tsne3d@cell.embeddings),file="tsne3d.txt",sep="\t",quote=F,row.names=F)

#MNN标化后，进行聚类分析

methodValue = "GraphCluster"
seuset <- FindNeighbors(seuset, reduction = "mnn", dims = pcadim)
seuset <- FindClusters(seuset, resolution = resolution)

table = table(Idents(seuset),seuset@meta.data$orig.ident)
table = dcast(data.frame(table),Var1~Var2)
colnames(table)[1]="Cluster"
write.table(table,file=paste0(methodValue,".Statistics.txt"),sep="\t",row.names=F,quot=F)
data = data.frame(Cell = colnames(seuset),Cluster=Idents(seuset))
write.table(data,file=paste0(methodValue,".Summary_Cell.txt"),sep="\t",row.names=F,quot=F)

#输出marker gene
allmarkers <- FindAllMarkers(object = seuset, only.pos = onlypos, min.pct = minpct,logfc.threshold = threshold, test.use = MarkerGeneMethod,base=base)
write.table(data.frame(allmarkers$gene,allmarkers[,1:6]),file=paste0(methodValue,".AllMarkerGenes.txt"),sep="\t",row.names=F,quot=F)

#保存seurat对象到rds
saveRDS(seuset, file = paste0("/MNN.seuset.rds"))


