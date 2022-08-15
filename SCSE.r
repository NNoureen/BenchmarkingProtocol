stringsAsFactors=FALSE
library(stringr)
library(foreach)
library(doParallel)
library('GSA')
registerDoParallel(20)  # use multicore, set to the number of our cores

################ SCSE function ############ 
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

##############################################################################################################################
############## Reading data into R  ############################################

data = readRDS('IDHAstrocytoma_GE_20210311.RDS')
data <- as.matrix(data)

############################# Reading genes sets and computing scores using SCSE function ################################

Genesets1 <-  GSA.read.gmt('h.all.v7.2.symbols.gmt')              ### loading hallmark gene sets
GSsize = length(Genesets1$genesets)


CombineSCSEResult1 <-  foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("C5")
	print(k)
    Result_SCSE1 = Execute_SCSE(data,Genesets1,k)          ### computing scores for hallmark gene sets using SCSE
	Result_SCSE1
}

write.table(CombineSCSEResult1,'SCSE_IDHA_HallmarksResult.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c2.all.v7.2.symbols.gmt')          ### loading C2 gene sets
GSsize = length(Genesets1$genesets)


CombineSCSEResult2 <-  foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("C5")
	print(k)
    Result_SCSE2 = Execute_SCSE(data,Genesets1,k)        ### computing scores for C2 gene sets using SCSE
	Result_SCSE2
}

write.table(CombineSCSEResult2,'SCSE_IDHA_C2_GSets.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c3.all.v7.2.symbols.gmt') ### loading C3 gene sets
GSsize = length(Genesets1$genesets)

CombineSCSEResult3 <-  foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("C5")
	print(k)
    Result_SCSE3 = Execute_SCSE(data,Genesets1,k)        ### computing scores for C3 gene sets using SCSE
	Result_SCSE3
}

write.table(CombineSCSEResult3,'SCSE_IDHA_C3.txt',sep='\t',quote=FALSE, row.names=FALSE)

