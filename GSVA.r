library(GSA)
library(GSVA)
stringsAsFactors=FALSE
library(stringr)
library(foreach)
library(doParallel)
registerDoParallel(200)  # use multicore, set to the number of our cores
####################################

Execute_GSVA <- function(exprMatrix,Genesets1,k){

genes = unlist(Genesets1$genesets[k])
geneslist = list(genes)
pathwayName = Genesets1$geneset.names[k]
Result = gsva(exprMatrix, geneslist,method="gsva")  ###GSVA 
Result = data.frame(Result)
Result$Msigb = pathwayName
return(Result)
}

################### Reading Data into R #################

HNSC = readRDS('IDHAstrocytoma_GE_20210311.RDS')
exprMatrix <- as.matrix(HNSC)

############################# Reading genes sets and calling gsva function ################################

Genesets1 <-  GSA.read.gmt('h.all.v7.2.symbols.gmt')         ### loading hallmark gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult1 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("GSVA:C5")
	print(k)
	Result_GSVA1  = Execute_GSVA(exprMatrix,Genesets1,k)              ### computing scores for hallmark gene sets using gsva
	Result_GSVA1
}

write.table(Combine_GSVAResult1,'GSVA_IDHA_HallmarksResult.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c2.all.v7.2.symbols.gmt')      ### loading C2 gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult2 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("GSVA:C5")
	print(k)
	Result_GSVA2  = Execute_GSVA(exprMatrix,Genesets1,k)         ### computing scores for C2 gene sets using gsva
	Result_GSVA2
}

write.table(Combine_GSVAResult2,'GSVA_IDHA_C2_GSets.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c3.all.v7.2.symbols.gmt')      ### loading C3 gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult3 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("GSVA:C5")
	print(k)
	Result_GSVA3  = Execute_GSVA(exprMatrix,Genesets1,k)          ### computing scores for C3 gene sets using gsva
	Result_GSVA3 
}

write.table(Combine_GSVAResult3,'GSVA_IDHA_C3.txt',sep='\t',quote=FALSE, row.names=FALSE)

