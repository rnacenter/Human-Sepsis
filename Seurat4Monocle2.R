rdsFile = "SeuratObject.rds"            #输入SeuratObject的rds
outpath = "Output_dir"                  #输出路径
geneFile = "genelist.txt"               #输入计算的基因，不输入使用variableGene
lambda = 2                              #DDRTree参数lambda，值越大节点越少

library(monocle)
library(Seurat)
library(reshape2)
library(dplyr)
setwd(outpath)
seuset = readRDS(rdsFile)

seuCDSAll = importCDS(seuset, import_all = TRUE)
if(!is.null(geneFile)){
	genelist = unique(read.table(geneFile,sep="\t",header=T)[,1])
}else{
	genelist = VariableFeatures(seuset)
}

seuCDSAll <- setOrderingFilter(seuCDSAll, genelist)
seuCDSAll <- estimateSizeFactors(seuCDSAll)
seuCDS <- reduceDimension(seuCDSAll, max_components = 2, method = "DDRTree",lambda = lambda*ncol(seuCDSAll))
seuCDS <- orderCells(seuCDS)

plot_cell_trajectory(seuCDS, color_by = "State")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
plot_cell_trajectory(seuCDS, color_by = "Cluster")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
plot_cell_trajectory(seuCDS, color_by = "Pseudotime")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

branchPoint = 1
numClusters = 3

BEAM_res <- BEAM(seuCDS, branch_point = branchPoint)
plot_genes_branched_heatmap(seuCDS,branch_point = branchPoint,num_clusters = numClusters,use_gene_short_name = T,show_rownames = T,return_heatmap=TRUE)
