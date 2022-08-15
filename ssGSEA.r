library('GSVA')
library('GSA')
stringsAsFactors=FALSE
library(circlize)
library(stringr)
library(foreach)
library(doParallel)
registerDoParallel(20)  # use multicore, set to the number of our cores
#################################### ssGSEA function 

Execute_GSVA <- function(exprMatrix,Genesets1,k){

genes = unlist(Genesets1$genesets[k])
geneslist = list(genes)
pathwayName = Genesets1$geneset.names[k]
Result = gsva(exprMatrix, geneslist,method="ssgsea")
Result = data.frame(Result)
Result$Msigb = pathwayName
return(Result)
}

#################################### Reading data into R

data = readRDS('IDHAstrocytoma_GE_20210311.RDS')
exprMatrix <- as.matrix(data)

############################# Reading genes sets and computing scores using ssGSEA ################################


Genesets1 <-  GSA.read.gmt('h.all.v7.2.symbols.gmt')            ### loading hallmark gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult1 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("GSVA:C5")
	print(k)
	Result_GSVA1  = Execute_GSVA(exprMatrix,Genesets1,k)        ### computing scores for hallmark gene sets using ssGSEA
	Result_GSVA1
}

write.table(Combine_GSVAResult1,'ssGSEA_IDHA_HallmarksResult.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c2.all.v7.2.symbols.gmt')        ### loading C2 gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult2 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("GSVA:C5")
	print(k)
	Result_GSVA2  = Execute_GSVA(exprMatrix,Genesets1,k)      ### computing scores for C2 gene sets using ssGSEA
	Result_GSVA2
}

write.table(Combine_GSVAResult2,'ssGSEA_IDHA_C2_GSets.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c3.all.v7.2.symbols.gmt')       ### loading C3 gene sets
GSsize = length(Genesets1$genesets)


Combine_GSVAResult3 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar% 
{
	print("GSVA:C5")
	print(k)
	Result_GSVA3  = Execute_GSVA(exprMatrix,Genesets1,k)       ### computing scores for C3 gene sets using ssGSEA
	Result_GSVA3
}

write.table(Combine_GSVAResult3,'ssGSEA_IDHA_C3.txt',sep='\t',quote=FALSE, row.names=FALSE)

