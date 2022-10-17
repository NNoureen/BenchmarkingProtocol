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


################################################ Signature Scoring of Random Up and DN regulated (GOLD STANDARD) gene sets ################



## Loading required package: iterators
registerDoParallel(20)  # use multicore, set to the number of our cores

Genesets1 <-  GSA.read.gmt('IDHAstrocytoma_RS_Up.gmt')                       ### reading Random gene sets for Up regulated genes
GSsize = length(Genesets1$genesets)


CombineJAS1 <- foreach (k=1:5000,.combine = rbind,.errorhandling = "remove") %dopar%   ## parallal processing using dopar
{
	print("C5")
	print(k)
	Result_JAS1 =     ExecutJAS(data,Genesets1,k)
	Result_JAS1
}

write.table(CombineJAS1,'JAS_Astrocytoma_UpRGS.txt',sep='\t',quote=FALSE, row.names=FALSE)   ### Saving the scores for Random Up regulated gene sets


Genesets1 <-  GSA.read.gmt('IDHAstrocytoma_DN.gmt')                    ### reading Random gene sets for DN regulated genes
GSsize = length(Genesets1$genesets)


CombineJAS2 <- foreach (k=1:5000,.combine = rbind,.errorhandling = "remove") %dopar%  ## parallal processing using dopar
{
	print("C5")
	print(k)
	Result_JAS1 =     ExecutJAS(data,Genesets1,k)
	Result_JAS1
}

write.table(CombineJAS2,'JAS_Astrocytoma_DNRGS.txt',sep='\t',quote=FALSE, row.names=FALSE)   ### Saving the scores for Random DN regulated gene sets
