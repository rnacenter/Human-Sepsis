library(GSVA)
library(limma)
library(stringr)

gmtFile = "./gene.gmt"
inputFile = "./expData.txt"

gmtread = function (file)
{
     if (!grepl("\\.gmt$", file)[1]) {
        data = read.table(file,header = T,sep="\t")
	types = as.character(unique(data[,1]))
	geneSetDB = list()
	for(i in 1:length(types)){
		type = types[i]
		geneset = as.character(unique(data[data[,1]==type,2]))
		geneSetDB[[i]]=geneset
	}
	names(geneSetDB) = types 
     }else{
     	geneSetDB = readLines(file)
     	geneSetDB = strsplit(geneSetDB, "\t")
     	names(geneSetDB) = sapply(geneSetDB, "[", 1)
    	geneSetDB = lapply(geneSetDB, "[", -1:-2)
     	geneSetDB = lapply(geneSetDB, function(x) {
        	x[which(x != "")]
     	})
     }
     return(geneSetDB)
}

geneset = gmtread(gmtFile)

data <- read.table(inputFile,sep="\t",header=T,check.names=F)

gsva_matrix<- gsva(data, geneset,method="ssgsea",ssgsea.norm=FALSE)
gsva_matrix = t(gsva_matrix)

write.table(gsva_matrix,file="ssgsea.txt",sep="\t",quot=F)

