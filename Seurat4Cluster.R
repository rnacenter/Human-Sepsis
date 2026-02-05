rdsFile = "Seurat_Norm_Object.rds"	#运行过Norm的SeuratObjectRDS
prefix = "object_name"				#输出名称前缀
outpath = "output_dir"				#输出结果路径
resolution = 0.8					#聚类分辨率

#计算FindAllMarkers参数
MarkerGeneMethod = "wilcox"     #计算markergene的算法
threshold = 0.25                #差异倍数 0.25以下不算显著Marker，会被去除，默认0.25
minpct = 0.1                    #计算MarkerGene表达占总体细胞的最低百分比
onlypos = TRUE                  #只返回positive基因
base=2                          #计算markergene的差异倍数使用2为底


library(plyr)
library(reshape2)
library(Seurat)
library(mclust)
library(dplyr)

setwd(outpath)
seuset = readRDS(rdsFile)
assay = DefaultAssay(seuset)

dir.create("GraphClust")
if(!is.null(seuset@commands$RunUMAP)){
	dims = seuset@commands$RunUMAP@params$dims
}else if(!is.null(seuset@commands$RunTSNE)){
	dims = seuset@commands$RunTSNE@params$dims
}
seuset <- FindNeighbors(seuset, reduction = "pca",dims = dims)
seuset <- FindClusters(seuset, resolution = resolution)
seuset <- BuildClusterTree(object = seuset)

table = table(Idents(seuset),seuset@meta.data$orig.ident)
table = dcast(data.frame(table),Var1~Var2)
colnames(table)[1]="Cluster"
write.table(table,file=paste0(prefix,".GraphClust.Statistics.txt"),sep="\t",row.names=F,quote=F)
data = data.frame(Cell = colnames(seuset),Cluster=Idents(seuset))
write.table(data,file=paste0(prefix,".GraphClust.Summary_Cell.txt"),sep="\t",row.names=F,quot=F)

allmarkers <- FindAllMarkers(object = seuset, only.pos = onlypos, min.pct = minpct,logfc.threshold = threshold, test.use = MarkerGeneMethod,base=base)
write.table(data.frame(allmarkers$gene,allmarkers[,1:6]),file=paste0(prefix,".GraphClust.AllMarkerGenes.txt"),sep="\t",row.names=F,quote=F)

DefaultAssay(seuset) <- assay
saveRDS(seuset, file = paste0(prefix,".GraphClust.seuset.rds"))

