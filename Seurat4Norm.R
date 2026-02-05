#利用Seurat对cellranger矩阵结果进行分析标化
#参数列表
countsFile = "input_count_dir"		#输入文件CellRanger的matrix结果，counts表的tsv或csv
genecolumn = 2						#选择matrix结果中 features表读第一列还是第二列
prefix = "object_name"				#输出名称前缀
outpath = "output_dir"				#输出结果路径
samplelist = NULL					#输入细胞样本信息(.txt)
celllist = NULL						#筛选输入文件中的细胞列表(.txt)
varGenelist = NULL					#根据输入文件指定varGene(.txt)
geneExpMinCell  = 3					#表示一个基因至少要在3个细胞里有表达量，否则过滤掉这个基因。
cellMinGene = 200					#每个细胞最少有200个基因表达，否则过滤掉这个细胞。
cellMaxGene = 10000					#表示一个细胞表达基因超过10000有可能是双细胞，被过滤掉
MaxMTPercent = 0.2					#表示一个细胞线粒体基因超过20%有可能是一个死细胞，会被过滤掉
MaxRBCPercent = 0.2					#表示一个细胞红细胞基因超过20%，会被过滤掉
MaxHouseKeepingPercent = 0.4		#表示一个细胞housekeeping基因超过40%，会被过滤掉
MaxDAPercent = 0.2					#表示一个细胞DissociationAssociated基因超过20%，会被过滤掉	
removeGene = "MTGene,HouseKeeping,DAgene,RBCGene,RPGene"		#选择要删除的基因
ExtraAnalysis = "MT,CellCycle,HouseKeeping,dagene,RBC,RPGene"	#选择要进行额外计算的基因
ScaleModel = 'linear'				#数据中心化处理的算法 Seurat提供的3种不同算法 默认linear (运算速度最快) poisson,negbinom;
NormType = 	'scale'					#标化方法： scale 或者 sctransform	
regressOut = "nCount_RNA,percent.mito"	#做线性回归时的参数设置	
#数据库路径(按需配置)：包含MT.txt,RBCGene.txt,RPGene.txt,HouseKeeping.txt,DAGene.txt,sGene.txt,g2mGene.txt
databasepath = "database_dir"

pcadim = 1:10					#手动输入的pca维度
nneighbors = 30					#UMAP分析neighbors数量
computePCs = 50					#主成分分析返回个数
varGeneNum = 2000					#手动输入的variablegene的数量

doubletfinder = TRUE				#是否载入DoubletFinder分析结果
doubletfinderRes = "DFRes.txt"		#载入doubletfinder结果
deleteDF = TRUE						#是否删除双细胞

library(plyr)
library(reshape2)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
options(stringsAsFactors = FALSE)

#默认3，表示一个基因至少要在三个细胞里有表达量，否则过滤掉这个基因。
geneExpMinCell = geneExpMinCell;
#默认200每个细胞最少有200个基因表达，否则过滤掉这个细胞。
cellMinGene = cellMinGene
#默认10000表示一个细胞表达基因超过5000有可能是双细胞，被过滤掉
cellMaxGene = cellMaxGene
#默认0.2表示线粒体基因超过20%有可能是一个死细胞，会被过滤掉
MaxMTPercent = MaxMTPercent
#默认0.2表示红细胞基因超过20%，会被过滤掉
MaxRBCPercent = MaxRBCPercent
#默认0.4 表示housekeeping基因超过40%，会被过滤掉
HouseKeepingPercent = MaxHouseKeepingPercent;
#默认0.2表示DissociationAssociated基因超过20%，会被过滤掉
maxDAPercent = MaxDAPercent;
#数据中心化处理的算法 Seurat提供的3种不同算法 默认linear 运算速度最快 poisson negbinom
ScaleModel = ScaleModel
#要删除的gene:MTGene;HouseKeeping Gene;Dissociation-Associated Gene;RBCGene
removeGene = unlist(strsplit(removeGene,"[,]"))
#选择要进行的额外分析：MT;CellCycle;HouseKeeping;Dissociation-Associated Gene;RBC
ExtraAnalysis = unlist(strsplit(ExtraAnalysis,"[,]"))
#手动输入的pca维度
pcadim = pcadim
#主成分分析返回个数
computePCs = computePCs
#手动输入的variablegene的数量
varGeneNum = varGeneNum
#UMAP分析neighbors数量
nneighbors = nneighbors
#标化方法： norm&Scale 或者 sctransform
NormType = NormType
#做线性回归时的参数设置：MTPerecnt,Sample,CountsNum,GeneNum,S.Score(CellCycleOnly)),G2M.Score(CellCycleOnly))
regressOut = regressOut
#使用默认参数还是手动设置参数
parameter = parameter
#对应物种数据库路径
databasepath = databasepath
#筛选输入文件中的细胞列表
celllist = celllist
#根据输入文件指定varGene
varGenelist = varGenelist
#选择matrix结果中 gene表读第一列还是第二列
genecolumn = genecolumn
doubletfinder = doubletfinder
#是否删除双细胞
deleteDF = deleteDF


setwd(outpath)

if(dir.exists(countsFile)){
    allcounts = Read10X(data.dir=countsFile,gene.column = genecolumn)
    if(class(allcounts)=="list"){
      	allcounts = allcounts[[1]]
    }
}
#如果输入细胞list，一列：取两者交集细胞，两列：取交集并添加初始ident数据
realOrigIdents = NULL
if(!is.null(celllist)){
    celllist = read.table(celllist,sep="\t",header=F,comment.char = "")
	rownames(celllist) = celllist[,1]
	realcells = intersect(colnames(allcounts),celllist[,1])
	allcounts = allcounts[,realcells]
	if(ncol(celllist)>1){
		realcells = celllist[realcells,]
		realOrigIdents = realcells[,2]
		names(realOrigIdents) = realcells[,1]
    }
}
#添加样本信息
if(!is.null(samplelist)){
    samplelist = samplelist;
    sampleInfo = read.table(samplelist, sep = "\t", header = FALSE,row.names=1)
    samples = sapply(strsplit(colnames(allcounts), "-"), function(v) return(v[2]))
    if(length(unique(samples))==1){
        samples[is.na(samples)]<-sampleInfo[1,1]
        sampleInfo = samples
    }else{
        sampleInfo = sampleInfo[samples,]
    }
    names(sampleInfo) = colnames(allcounts)
    #seuset <- AddMetaData(object = seuset, metadata = sampleInfo, col.name = "orig.ident")
}else{
	sampleInfo = rep(prefix,ncol(allcounts))
	names(sampleInfo) = colnames(allcounts)
}

dir.create("1.Normalization")
dir.create("2.PCAAnalysis")

#将counts表中表达量为0的列删掉
allcounts = allcounts[,colSums(allcounts)>0,drop=F]

#根据选择的额外分析分情况
if("MT" %in% ExtraAnalysis | "MTGene" %in% removeGene){
	mito.genes = read.table(paste0(databasepath,"/MT.txt"), sep = "\t", header = FALSE)[,1]
	mito.genes = mito.genes[is.finite(match(mito.genes,rownames(allcounts)))]
	print(length(mito.genes))
	percent.mito <- Matrix::colSums(allcounts[mito.genes,,drop=F])/Matrix::colSums(allcounts)
}
if("RBC" %in% ExtraAnalysis | "RBCGene" %in% removeGene){
	rbc.genes = read.table(paste0(databasepath,"/RBCGene.txt"), sep = "\t", header = FALSE)[,1]
	rbc.genes = rbc.genes[is.finite(match(rbc.genes,rownames(allcounts)))]
	print(length(rbc.genes))
	percent.rbc <- Matrix::colSums(allcounts[rbc.genes,,drop=F])/Matrix::colSums(allcounts)

}
if("RPGene" %in% ExtraAnalysis | "RPGene" %in% removeGene){
	rp.genes = read.table(paste0(databasepath,"/RPGene.txt"), sep = "\t", header = FALSE)[,1]
	rp.genes = rp.genes[is.finite(match(rp.genes,rownames(allcounts)))]
	print(length(rp.genes))
	percent.rp <- Matrix::colSums(allcounts[rp.genes,,drop=F])/Matrix::colSums(allcounts)
}

if("HouseKeeping" %in% ExtraAnalysis | "HouseKeeping" %in% removeGene){
    HouseKeepinggenes = read.table(paste0(databasepath,"/HouseKeeping.txt"), sep = "\t", header = FALSE)[,1]
	num = length(HouseKeepinggenes)
	HouseKeepinggenes = HouseKeepinggenes[is.finite(match(HouseKeepinggenes,rownames(allcounts)))]
	print(length(HouseKeepinggenes))
	percent.HouseKeeping = apply(allcounts[HouseKeepinggenes, ],2,FUN=function(x){return(sum(x > 0) / num)})
}

if("dagene" %in% ExtraAnalysis | "DAgene" %in% removeGene){

	dagenes = read.table(paste0(databasepath,"/DAGene.txt"), sep = "\t", header = FALSE)[,1]
	dagenes = dagenes[is.finite(match(dagenes,rownames(allcounts)))]
	print(length(dagenes))
	percent.dagene <- Matrix::colSums(allcounts[dagenes,,drop=F])/Matrix::colSums(allcounts)
}


if(doubletfinder){
	doubletresult <-read.table(doubletfinderRes,sep="\t",header=T,check.names=F)
}
#构建seuratobject
seuset <- CreateSeuratObject(counts = allcounts,project = prefix,min.cells = geneExpMinCell, min.features = 1)
if("MTGene" %in% removeGene){
    print("Remove MTGene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),mito.genes))
}

if("RBCGene" %in% removeGene){
    print("Remove RBCGene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),rbc.genes))
}
if("RPGene" %in% removeGene){
    print("Remove RPGene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),rp.genes))
}
if("HouseKeeping" %in% removeGene){
    print("Remove HouseKeeping Gene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),HouseKeepinggenes))
}
if("DAgene" %in% removeGene){
    print("Remove DissociationAssociated Gene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),dagenes))
}
print("nFeature_RNA")
quantile(seuset@meta.data$nFeature_RNA,probs = seq(0, 1, 0.05))
if("MT" %in% ExtraAnalysis){
	seuset <- AddMetaData(object = seuset, metadata = percent.mito, col.name = "percent.mito")
}else{
	seuset <- AddMetaData(object = seuset, metadata = 0, col.name = "percent.mito")
}

if("RBC" %in% ExtraAnalysis){
	seuset <- AddMetaData(object = seuset, metadata = percent.rbc, col.name = "percent.rbc")
}

if("RPGene" %in% ExtraAnalysis){
	seuset <- AddMetaData(object = seuset, metadata = percent.rp, col.name = "percent.rp")
}
#添加样本信息
if(!is.null(sampleInfo)){
    seuset <- AddMetaData(object = seuset, metadata = sampleInfo, col.name = "orig.ident")
}
#添加ident信息
if(!is.null(realOrigIdents)){
    seuset <- AddMetaData(object = seuset, metadata = realOrigIdents, col.name = "orig.ident")
}
nSample = length(unique(seuset@meta.data$orig.ident))

if("HouseKeeping" %in% ExtraAnalysis){
    seuset <- AddMetaData(object = seuset, metadata = percent.HouseKeeping, col.name = "percent.HouseKeeping")
}
if("dagene" %in% ExtraAnalysis){
    seuset <- AddMetaData(object = seuset, metadata = percent.dagene, col.name = "percent.DissociationAssociated")
}
if(doubletfinder){
	singletinfo <- doubletresult[is.finite(match(doubletresult[,1],colnames(seuset))),,drop=F]
	seuset <- AddMetaData(object = seuset, metadata = as.numeric(singletinfo[,2]), col.name = "DoubletValue")
	seuset <- AddMetaData(object = seuset, metadata = singletinfo[,3], col.name = "DoubletClassification")
}
write.table(data.frame(CellName=rownames(seuset@meta.data),seuset@meta.data),file="1.Normalization/CellInfosRaw.txt",sep="\t",quote=F,row.names=F)

#删除双细胞
if(doubletfinder){
	if(deleteDF){
		singlet <- doubletresult[which(doubletresult[,3]=="Singlet"),1]
		seuset <- subset(seuset,cells=as.character(singlet))
	}
}
#添加新的ident
seuset = SetIdent(object = seuset, cells = colnames(seuset), value = seuset@meta.data$orig.ident)
#筛选满足nGene区间的细胞
print(paste0("nGeneFilterRegion:",cellMinGene,"-",cellMaxGene))
seuset <- subset(x = seuset, subset = nFeature_RNA > cellMinGene & nFeature_RNA < cellMaxGene)
#细胞条件筛选
if("MT" %in% ExtraAnalysis){
	print(paste0("MTPercentFilterRegion:<",MaxMTPercent))
	seuset <- subset(x = seuset, subset = percent.mito < MaxMTPercent)
}
if("RBC" %in% ExtraAnalysis){
	print(paste0("RBCPercentFilterRegion:<",MaxRBCPercent))
	seuset <- subset(x = seuset, subset = percent.rbc < MaxRBCPercent)
}

if("HouseKeeping" %in% ExtraAnalysis){
    print(paste0("HouseKeepingPercentFilterRegion:>=",HouseKeepingPercent))
	seuset <- subset(x = seuset, subset = percent.HouseKeeping >= HouseKeepingPercent) 
}
if("dagene" %in% ExtraAnalysis){
    print(paste0("Dissociation-AssociatedPercentFilterRegion:<",maxDAPercent))
	seuset <- subset(x = seuset, subset = percent.DissociationAssociated < maxDAPercent)
}


seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", scale.factor = 10000)

seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = varGeneNum)

if("CellCycle" %in% ExtraAnalysis){
	Sgene <-read.table(paste0(databasepath,"/sGene.txt"), sep = "\t", header = FALSE)[,1]
	G2Mgene <- read.table(paste0(databasepath,"/g2mGene.txt"), sep = "\t", header = FALSE)[,1]
	seuset = CellCycleScoring(seuset, s.features =Sgene, g2m.features = G2Mgene,  set.ident = FALSE)
	seuset$CC.Difference <- seuset$S.Score - seuset$G2M.Score
}
if(("dagene" %in% ExtraAnalysis)&& (!"DAgene" %in% removeGene)){
	seuset = AddModuleScore(seuset,list(dagenes),name="DissociationAssociated.Score")
	colnames(seuset@meta.data)[colnames(seuset@meta.data)=="DissociationAssociated.Score1"]="DissociationAssociated.Score"
}

if(!is.null(regressOut)){
    regressOut = unlist(strsplit(regressOut, "[,]"))
    if(length(unique(seuset@meta.data$orig.ident))==1){
        regressOut = regressOut[!regressOut%in%"orig.ident"]
    }
    if(!"CellCycle" %in% ExtraAnalysis){
        regressOut = regressOut[!regressOut%in%c("S.Score","G2M.Score","CC.Difference")]
    }
	if(!"dagene" %in% ExtraAnalysis){
        regressOut = regressOut[!regressOut%in%c("percent.DissociationAssociated","DissociationAssociated.Score")]
    }
    if(NormType=="sctransform"){
        regressOut = regressOut[!regressOut%in%c("nFeature_RNA","nCount_RNA")]
    }
    print(regressOut)
}

#如果输入varGene列表则将输入基因作为varGene
if(!is.null(varGenelist)){
    varGenelist = read.table(varGenelist,sep="\t",header=F,comment.char = "")[,1]
    varGenelist = intersect(rownames(seuset),varGenelist)
    VariableFeatures(seuset) = varGenelist
}

#标化为norm&scale时，scaledata
if(NormType=="scale"){
	write.table(VariableFeatures(seuset),file="1.Normalization/varGene.txt",sep="\t",row.names=F,quote=F)
	if(RunFast){
        scalegene = VariableFeatures(seuset)
        print(paste0("Use Variable genes to scale data"))
	}else{
        scalegene = rownames(seuset)
	}
	seuset <- ScaleData(object = seuset, features = scalegene, vars.to.regress = regressOut, model.use = ScaleModel)
}
#标化类型为sctransform，运行SCTransform，assay为SCT
if(NormType=="sctransform"){
	seuset <- SCTransform(seuset, vars.to.regress = regressOut, variable.features.n = varGeneNum, return.only.var.genes=RunFast)
	write.table(VariableFeatures(seuset),file="1.Normalization/varGene.txt",sep="\t",row.names=F,quote=F)
}
write.table(data.frame(CellName=rownames(seuset@meta.data),seuset@meta.data),file="1.Normalization/CellInfosFilter.txt",sep="\t",quote=F,row.names=F)

seuset <- RunPCA(object = seuset, features = VariableFeatures(seuset), npcs = computePCs, do.print = TRUE, ndims.print = 1:computePCs, nfeatures.print = 10,seed.use=2026)
if(NormType=="scale"){
	num.replicate = 100
	seuset <- JackStraw(seuset, num.replicate = num.replicate, dims = computePCs)
	seuset <- ScoreJackStraw(seuset, dims = 1:computePCs)
	print(seuset@reductions$pca@jackstraw@overall.p.values)
	write.table(seuset@reductions$pca@jackstraw@overall.p.values,file="2.PCAAnalysis/pcapvalue.txt",sep="\t",quote=F,row.names=F)
}
#计算合并pcgene
pcgene <- c()
for(i in 1:computePCs){
    topgenes = TopFeatures(object = seuset[["pca"]], dim = i,balanced=TRUE, nfeatures=60)
    topgenes = data.frame(PC=paste0("PC",i),topgenes$positive,topgenes$negative)
    pcgene = rbind(pcgene,topgenes)
}
write.table(pcgene,file="2.PCAAnalysis/pcGene.txt",sep="\t",quote=F,row.names=F,col.names=T)
write.table(data.frame(CellName=rownames(seuset@reductions$pca@cell.embeddings),seuset@reductions$pca@cell.embeddings),file="2.PCAAnalysis/pca.txt",sep="\t",quote=F,row.names=F)

print("RunUMAP")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors, reduction.name = "umap3d", n.components = 3, reduction.key = "umap3d_")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors)
#写出umap坐标
write.table(data.frame(CellName=rownames(seuset@reductions$umap@cell.embeddings),seuset@reductions$umap@cell.embeddings),file=paste0(prefix,".umapCoords.2d.txt"),sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$umap3d@cell.embeddings),seuset@reductions$umap3d@cell.embeddings),file=paste0(prefix,".umapCoords.3d.txt"),sep="\t",quote=F,row.names=F)

print("RunTSNE")
seuset <- RunTSNE(object = seuset,dims = pcadim,reduction.name = "tsne3d",dim.embed=3, reduction.key = "tSNE3d_",check_duplicates = FALSE)
seuset <- RunTSNE(object = seuset,dims = pcadim,check_duplicates = FALSE)
#写出tsne坐标
write.table(data.frame(CellName=rownames(seuset@reductions$tsne@cell.embeddings),seuset@reductions$tsne@cell.embeddings),file=paste0(prefix,".tsneCoords.2d.txt"),sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$tsne3d@cell.embeddings),seuset@reductions$tsne3d@cell.embeddings),file=paste0(prefix,".tsneCoords.3d.txt"),sep="\t",quote=F,row.names=F)

#保存RDS	
print("SaveRDS")
saveRDS(seuset, file = paste0(prefix,".SeuratObject.rds"))


