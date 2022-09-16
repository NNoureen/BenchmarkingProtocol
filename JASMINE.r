stringsAsFactors=FALSE
library(stringr)
library(foreach)
library(doParallel)
library('GSA')          ## for calling GSA.read.gmt function
registerDoParallel(20)  # use multicore, set to the number of cores to 20


################################## JASMINE function #########################################################

### Rank Calculation function in JASMINE
## This function is called by JASMINE function to compute ranks of each cell for a gene signature
 
RankCalculation <- function(x,genes){
			##calculate relative rank of expressed genes from the gene signature
            subdata = x[x!=0]
			DataRanksUpdated=rank(subdata)
			DataRanksSigGenes = DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)]
			CumSum = ifelse(length(DataRanksSigGenes),mean(DataRanksSigGenes,na.rm = TRUE),0 ) ### Cumsum of each sample calculated based on the ranks
			FinalRawRank = CumSum/length(subdata)
			return(FinalRawRank)
			}			

### OR Calculation function in JASMINE
## This function is called by JASMINE function to compute enrichment of a gene signature per cell
 
			
ORCalculation <- function(data,genes){
			GE = data[which(rownames(data) %in% genes),]
			NGE = data[-which(rownames(data) %in% genes),]
			SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))
			NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))
			SigGenesNE = nrow(GE) - SigGenesExp
			SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)			
			NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)
		    NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)
			NSigGenesNE = NSigGenesNE - SigGenesNE
			OR = (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp)
            return(OR)
			}

### Normalization function in JASMINE
## This function is called by JASMINE function to scale the scores from 0 to 1 for a gene set across all the cells
 
NormalizationJAS <- function(JAS_Scores)
            {
				JAS_Scores = (JAS_Scores - min(JAS_Scores))/(max(JAS_Scores)- min(JAS_Scores))
				return(JAS_Scores)
			}


## JASMINE function takes as input the single-cell gene expression data, genes in a gene set and pathwayName/gene set name
 			
			
JASMINE <- function(data,genes,pathwayName)
		{
  		    idx = match(genes,rownames(data))
	        idx = idx[!is.na(idx)]
			if(length(idx)> 1){
			RM = apply(data,2,function(x) RankCalculation(x,genes))
			OR = ORCalculation(data,genes)			
			RM = NormalizationJAS(RM)
            OR = NormalizationJAS(OR)
			JAS_Scores = (RM + OR)/2
			FinalScores = data.frame(names(RM),JAS_Scores)
			colnames(FinalScores)[1]='SampleID'
			Samples = FinalScores$SampleID
			FinalScores = FinalScores[,-1]
			FinalScores = data.frame(t(FinalScores))
			names(FinalScores)=as.character(Samples)
			FinalScores$Msigb = pathwayName
			
			return(FinalScores)
			}
		}

############################# Execution of JASMINE function

## Function Input:  single-cell gene expression data, all gene sets in a list, and k represents the index for a gene set in a list

## Function Output: Scores for each cell/sample ID along with pathwayName

ExecutJAS <- function(data,Genesets1,k){

genes = unlist(Genesets1$genesets[k])
pathwayName = Genesets1$geneset.names[k]
subsetJAS = JASMINE(data,genes,pathwayName)
return(subsetJAS)
}


###################### Read data #####################

data = readRDS('IDHAstrocytoma_GE_20210311.RDS')   ## Reading gene expression data in RDS format
data <- as.matrix(data)

################################################ Execute JASMINE using Msigb gene sets and dopar function ################


Genesets1 <-  GSA.read.gmt('h.all.v7.2.symbols.gmt')       ### Reading hallmark gene sets
GSsize = length(Genesets1$genesets)

CombineAUCell1 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%           ## parallal processing using dopar
{
	print("C5")
	print(k)
	Result_AUCell1 =    ExecutJAS(data,Genesets1,k)  ### computing scores for hallmark gene sets using JASMINE
	Result_AUCell1
}

write.table(CombineAUCell1,'JAS_IDHA_HallmarksResult.txt',sep='\t',quote=FALSE, row.names=FALSE)


Genesets1 <-  GSA.read.gmt('c2.all.v7.2.symbols.gmt')        ### Reading C2 gene sets
GSsize = length(Genesets1$genesets)

CombineAUCell2 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%               ## parallal processing using dopar
{
	print("C5")
	print(k)
	Result_AUCell2 =    ExecutJAS(data,Genesets1,k)                 ### computing scores for C2 gene sets using JASMINE
	Result_AUCell2
}

write.table(CombineAUCell2,'JAS_IDHA_C2_GSets.txt',sep='\t',quote=FALSE, row.names=FALSE)

Genesets1 <-  GSA.read.gmt('c3.all.v7.2.symbols.gmt')                ### Reading C3 gene sets
GSsize = length(Genesets1$genesets)

CombineAUCell3 <- foreach (k=1:GSsize,.combine = rbind,.errorhandling = "remove") %dopar%               ## parallal processing using dopar
{
	print("C5")
	print(k)
	Result_AUCell3 =    ExecutJAS(data,Genesets1,k)                    ### computing scores for C3 gene sets using JASMINE
	Result_AUCell3
}

write.table(CombineAUCell3,'JAS_IDHA_C3_GSets.txt',sep='\t',quote=FALSE, row.names=FALSE)

