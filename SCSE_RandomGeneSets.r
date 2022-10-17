stringsAsFactors=FALSE
library(stringr)
library(foreach)
library(doParallel)
library('GSA')         ## for calling GSA.read.gmt function
registerDoParallel(20)  # use multicore, set the number of cores to 20

################ SCSE function ############ 

## Function Input: single-cell gene expression data and genes in one gene sets (This function is implemented using the equation described in Pont et al.,2019)

## Function Output: Scores for each cell/sample ID


SingleCellSigExplorer <- function(data,genes)
{
	DataRanks = data[which(rownames(data) %in% genes),]
	if(length(nrow(DataRanks))!=0)
	{
		CumSum = data.frame(colSums(DataRanks, na.rm = FALSE, dims = 1))
		colnames(CumSum)[1]='RawRankSum'
		SampleID = rownames(CumSum)
		CumSum1 = data.frame(SampleID,CumSum)
		row.names(CumSum1)=NULL
		TotalUMICount = data.frame(colSums(data, na.rm = FALSE, dims = 1))
		colnames(TotalUMICount)[1]='TotalUMISum'
		SampleID = rownames(TotalUMICount)
		TotalUMICount = data.frame(SampleID,TotalUMICount)
		row.names(TotalUMICount)=NULL
		FinalScore = merge(CumSum1, TotalUMICount, by='SampleID')
		FinalScore$SingleScorer = FinalScore$RawRankSum/FinalScore$TotalUMISum*100
		colnames(FinalScore)[4] = 'scSigExp'
		FinalScore = FinalScore[,c(1,4)]
		return(FinalScore)
	} 
}

############### Calling SCSE function using gene sets #########################

## Function Input: single-cell gene expression data, all gene sets in a list, and k represents the index for a gene set in a list

## Function Output: A dataframe containing sample/cell ID, gene set score and pathwayName



Execute_SCSE <- function(data,Genesets1,k)
{
	genes = unlist(Genesets1$genesets[k])
	pathwayName = Genesets1$geneset.names[k]
	SCSE = SingleCellSigExplorer(data,genes)
	if(length(SCSE)!=0)
	{
		Samples = SCSE$SampleID
		SCSE = data.frame(SCSE[,-1])
		SCSE = data.frame(t(SCSE))
		names(SCSE)= Samples
		SCSE$Msigb= pathwayName
		return(SCSE)
	}
}


################################################ Signature Scoring of Random Up and DN regulated (GOLD STANDARD) gene sets ################

data = readRDS('IDHAstrocytoma_GE_20210311.RDS')             ## Reading gene expression data in RDS format
data <- as.matrix(data)
#####################################################################################################
## Loading required package: iterators
registerDoParallel(20)  # use multicore, set to the number of our cores

Genesets1 <-  GSA.read.gmt('IDHAstrocytoma_RS_Up.gmt')                       ### reading Random gene sets for Up regulated genes
GSsize = length(Genesets1$genesets)


CombineSCSEResult1 <-  foreach (k=1:5000,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("C5")
	print(k)
    Result_SCSE1 = Execute_SCSE(data,Genesets1,k)
	Result_SCSE1
}

write.table(CombineSCSEResult1,'SCSE_Astrocytoma_UpRGS1.txt',sep='\t',quote=FALSE, row.names=FALSE)    ### Saving the scores for Random UP regulated gene sets


Genesets1 <-  GSA.read.gmt('IDHAstrocytoma_DN.gmt')                    ### reading Random gene sets for DN regulated genes
GSsize = length(Genesets1$genesets)


CombineSCSEResult2 <-  foreach (k=1:5000,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("C5")
	print(k)
    Result_SCSE2 = Execute_SCSE(data,Genesets1,k)
	Result_SCSE2
}

write.table(CombineSCSEResult2,'SCSE_Astrocytoma_DNRGS.txt',sep='\t',quote=FALSE, row.names=FALSE)   ### Saving the scores for Random DN regulated gene sets

