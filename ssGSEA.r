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

############################# Reading genes sets and computing scores using ssGSEA and dopar function ################################


Genesets1 <-  GSA.read.gmt('h.all.v7.2.symbols.gmt')            ### Reading hallmark gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult1 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%  ## parallal processing using dopar
{
	print("GSVA:C5")
	print(k)
	Result_GSVA1  = Execute_GSVA(exprMatrix,Genesets1,k)        ### computing scores for hallmark gene sets using ssGSEA
	Result_GSVA1
}

write.table(Combine_GSVAResult1,'ssGSEA_IDHA_HallmarksResult.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c2.all.v7.2.symbols.gmt')        ### Reading C2 gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult2 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%      ## parallal processing using dopar
{
	print("GSVA:C5")
	print(k)
	Result_GSVA2  = Execute_GSVA(exprMatrix,Genesets1,k)      ### computing scores for C2 gene sets using ssGSEA
	Result_GSVA2
}

write.table(Combine_GSVAResult2,'ssGSEA_IDHA_C2_GSets.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c3.all.v7.2.symbols.gmt')       ### Reading C3 gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult3 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%      ## parallal processing using dopar
{
	print("GSVA:C5")
	print(k)
	Result_GSVA3  = Execute_GSVA(exprMatrix,Genesets1,k)       ### computing scores for C3 gene sets using ssGSEA
	Result_GSVA3
}

write.table(Combine_GSVAResult3,'ssGSEA_IDHA_C3.txt',sep='\t',quote=FALSE, row.names=FALSE)

