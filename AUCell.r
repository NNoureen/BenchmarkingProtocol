stringsAsFactors=FALSE
library(stringr)
library(AUCell)          ## provide AUCell function
library(foreach)
library(doParallel)
library('GSA')          ## for calling GSA.read.gmt function
registerDoParallel(20)  # use multicore, set to the number of our cores
########################### AUCell Function ###############################################
AUCellfunc <- function(exprMatrix,genes){

out <- tryCatch(
{
cells_rankings <- AUCell_buildRankings(exprMatrix,plotStats=FALSE)
geneSets <- list(geneSet1=genes)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cellsAUCretrieve = getAUC(cells_AUC)
cells_AUCellScores = t(cellsAUCretrieve)
cells_AUCellScores = data.frame(rownames(cells_AUCellScores),cells_AUCellScores)
colnames(cells_AUCellScores) = c('SampleID','AUCell')
row.names(cells_AUCellScores)= NULL
return(cells_AUCellScores)

},

error = function(cond) {
cells_AUCellScores = data.frame(SampleID = "NA", AUCell = "NA")
return(cells_AUCellScores)
}
)
return(out)
}

##############################################################################################################################

Execute_AUCell <- function(exprMatrix,Genesets1,k)
{

	genes = unlist(Genesets1$genesets[k])
	pathwayName = Genesets1$geneset.names[k]
	AUCellMethod = AUCellfunc(exprMatrix,genes)
	if(nrow(AUCellMethod) == ncol(exprMatrix))
	{
		Samples = AUCellMethod$SampleID
		AUCellMethod = data.frame(AUCellMethod[,-1])
		AUCellMethod = data.frame(t(AUCellMethod))
		names(AUCellMethod) = Samples
		AUCellMethod$Msigb = pathwayName
		return(AUCellMethod)
	}
}

##############  Reading Data set ############################################
HNSC = readRDS('IDHAstrocytoma_GE_20210311.RDS')
exprMatrix <- as.matrix(HNSC)

############################# Reading genes sets and Calling AUCell function ################################


Genesets1 <-  GSA.read.gmt('h.all.v7.2.symbols.gmt')                ### loading hallmark gene sets
GSsize = length(Genesets1$genesets)

CombineAUCell1 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("C5")
	print(k)
	Result_AUCell1 =    Execute_AUCell(exprMatrix,Genesets1,k)        ### computing scores for hallmark gene sets using AUCell
	Result_AUCell1
}

write.table(CombineAUCell1,'AUCell_IDHA_HallmarksResult.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c2.all.v7.2.symbols.gmt')             ### loading C2 gene sets
GSsize = length(Genesets1$genesets)

CombineAUCell2 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("C5")
	print(k)
	Result_AUCell2 =    Execute_AUCell(exprMatrix,Genesets1,k)         ### computing scores for C2 gene sets using AUCell
	Result_AUCell2
}

write.table(CombineAUCell2,'AUCell_IDHA_C2_GSets.txt',sep='\t',quote=FALSE, row.names=FALSE)

Genesets1 <-  GSA.read.gmt('c3.all.v7.2.symbols.gmt')           ### loading C3 gene sets
GSsize = length(Genesets1$genesets)


CombineAUCell3 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("C5")
	print(k)
	Result_AUCell3 =    Execute_AUCell(exprMatrix,Genesets1,k)      ### computing scores for C3 gene sets using AUCell
	Result_AUCell3
}

write.table(CombineAUCell3,'AUCell_IDHA_C3.txt',sep='\t',quote=FALSE, row.names=FALSE)

