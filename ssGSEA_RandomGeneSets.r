library('GSVA')          ## for calling gsva function
library('GSA')           ## for calling GSA.read.gmt function
stringsAsFactors=FALSE
library(circlize)
library(stringr)
library(foreach)
library(doParallel)
registerDoParallel(20)  # use multicore, set the number of cores to 20
#################################### ssGSEA function 

## Function Input: single-cell gene expression data, all gene sets in a list, and k represents the index for a gene set in a list

## Function Output: Scores for each cell/sample ID along with pathwayName


Execute_GSVA <- function(exprMatrix,Genesets1,k){

genes = unlist(Genesets1$genesets[k])
geneslist = list(genes)
pathwayName = Genesets1$geneset.names[k]
Result = gsva(exprMatrix, geneslist,method="ssgsea")      ### GSVA Function calling specifying ssGSEA as a method
Result = data.frame(Result)
Result$Msigb = pathwayName
return(Result)
}

#################################### Reading data into R

data = readRDS('IDHAstrocytoma_GE_20210311.RDS')                 ## Reading gene expression data in RDS format
exprMatrix <- as.matrix(data)


################################################ Signature Scoring of Random Up and DN regulated (GOLD STANDARD) gene sets ################


## Loading required package: iterators
registerDoParallel(20)  # use multicore, set to the number of our cores

Genesets1 <-  GSA.read.gmt('IDHAstrocytoma_RS_Up.gmt')                       ### reading Random gene sets for Up regulated genes
GSsize = length(Genesets1$genesets)


Combine_GSVAResult1 <- foreach (k=1:5000,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("GSVA:C5")
	print(k)
	Result_GSVA1  = Execute_GSVA(exprMatrix,Genesets1,k)
	Result_GSVA1
}

write.table(Combine_GSVAResult1,'ssGSEA_Astrocytoma_UpRGS.txt',sep='\t',quote=FALSE, row.names=FALSE)   ### Saving the scores for Random UP regulated gene sets



Genesets1 <-  GSA.read.gmt('IDHAstrocytoma_DN.gmt')                    ### reading Random gene sets for DN regulated genes
GSsize = length(Genesets1$genesets)



Combine_GSVAResult2 <- foreach (k=1:5000,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("GSVA:C5")
	print(k)
	Result_GSVA1  = Execute_GSVA(exprMatrix,Genesets1,k)
	Result_GSVA1
}

write.table(Combine_GSVAResult2,'ssGSEA_Astrocytoma_DNRGS.txt',sep='\t',quote=FALSE, row.names=FALSE)   ### Saving the scores for Random DN regulated gene sets
