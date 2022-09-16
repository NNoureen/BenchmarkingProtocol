library(GSA)                 ## for calling GSA.read.gmt function
library(GSVA)
stringsAsFactors=FALSE
library(stringr)
library(foreach)
library(doParallel)          ## for parallal processing
registerDoParallel(200)  # use multicore, set to the number of 200 cores


#################################### GSVA function calling using GSVA package ##############################

## Function Input: single-cell gene expression data, all gene sets in a list, and k represents the index for a gene set in a list

## Function Output: Scores for each cell/sample ID along with pathwayName



Execute_GSVA <- function(exprMatrix,Genesets1,k){

genes = unlist(Genesets1$genesets[k])
geneslist = list(genes)
pathwayName = Genesets1$geneset.names[k]
Result = gsva(exprMatrix, geneslist,method="gsva")  ### GSVA Function calling specifying gsva as a method
Result = data.frame(Result)
Result$Msigb = pathwayName
return(Result)
}

################### Reading Data into R #################

HNSC = readRDS('IDHAstrocytoma_GE_20210311.RDS')            ## Reading gene expression data in RDS format
exprMatrix <- as.matrix(HNSC)

############################# Reading genes sets and calling gsva function ################################

Genesets1 <-  GSA.read.gmt('h.all.v7.2.symbols.gmt')         ### Reading hallmark gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult1 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%       ## parallal processing using dopar
{
	print("GSVA:C5")
	print(k)
	Result_GSVA1  = Execute_GSVA(exprMatrix,Genesets1,k)              ### computing scores for hallmark gene sets using gsva
	Result_GSVA1
}

write.table(Combine_GSVAResult1,'GSVA_IDHA_HallmarksResult.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c2.all.v7.2.symbols.gmt')      ### Reading C2 gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult2 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%              ## parallal processing using dopar
{
	print("GSVA:C5")
	print(k)
	Result_GSVA2  = Execute_GSVA(exprMatrix,Genesets1,k)         ### computing scores for C2 gene sets using gsva
	Result_GSVA2
}

write.table(Combine_GSVAResult2,'GSVA_IDHA_C2_GSets.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c3.all.v7.2.symbols.gmt')      ### Reading C3 gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult3 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%           ## parallal processing using dopar
{
	print("GSVA:C5")
	print(k)
	Result_GSVA3  = Execute_GSVA(exprMatrix,Genesets1,k)          ### computing scores for C3 gene sets using gsva
	Result_GSVA3 
}

write.table(Combine_GSVAResult3,'GSVA_IDHA_C3.txt',sep='\t',quote=FALSE, row.names=FALSE)

