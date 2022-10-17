stringsAsFactors=FALSE
library(stringr)
library(AUCell)          ## provide AUCell function
library(foreach)        
library(doParallel)     ## for parallal processing
library('GSA')          ## for calling GSA.read.gmt function
registerDoParallel(20)  # use multicore, set to the number of cores

########################### AUCell Function to compute the scores ###############################################

## Function Input: single-cell gene expression data and genes in one gene sets

## Function Output: Scores for each cell/sample ID

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

############################## Calling of AUCell function using the Random Gene Sets ####################################

## Function Input: single-cell gene expression data, all gene sets in a list, and k represents the index for a gene set in a list

## Function Output: A dataframe containing sample/cell ID, gene set score and pathwayName


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

##############  Reading the single-cell gene expression Data set ############################################

HNSC = readRDS('IDHAstrocytoma_GE_20210311.RDS')             ## Reading gene expression data in RDS format
exprMatrix <- as.matrix(HNSC)


######################################## Signature Scoring of Random Up and DN regulated (GOLD STANDARD) gene sets ###################################


## Loading required package: iterators
registerDoParallel(20)  # use multicore, set to the number of our cores

Genesets1 <-  GSA.read.gmt('IDHAstrocytoma_RS_Up.gmt')                       ### reading Random gene sets for Up regulated genes
GSsize = length(Genesets1$genesets)

CombineAUCell1 <- foreach (k=1:5000,.combine = rbind,.errorhandling = "remove") %dopar%   ## parallal processing using dopar
{
	print("C5")
	print(k)
	Result_AUCell1 =    Execute_AUCell(exprMatrix,Genesets1,k)
	Result_AUCell1
}

write.table(CombineAUCell1,'AUCell_Astrocytoma_UpRGS.txt',sep='\t',quote=FALSE, row.names=FALSE)   ### Saving the scores for Random Up regulated gene sets



Genesets1 <-  GSA.read.gmt('IDHAstrocytoma_DN.gmt')                    ### reading Random gene sets for DN regulated genes
GSsize = length(Genesets1$genesets)

CombineAUCell2 <- foreach (k=1:5000,.combine = rbind,.errorhandling = "remove") %dopar%   ## parallal processing using dopar
{
	print("C5")
	print(k)
	Result_AUCell1 =    Execute_AUCell(exprMatrix,Genesets1,k)
	Result_AUCell1
}

write.table(CombineAUCell2,'AUCell_Astrocytoma_DNRGS.txt',sep='\t',quote=FALSE, row.names=FALSE)  ### Saving the scores for Random DN regulated gene sets
