library('dplyr')
stringsAsFactors=FALSE
library(stringr)
library(data.table)
library(effectsize)      ## include effectsize library

############################## Effect Size Function Definition #########################################


Metadata = read.table('IDHAstrocytoma_Metadata.txt',sep='\t',head=T)              ## Reading the Metadata to divide cells into tumor and normal cases
TumorCells = Metadata$SampleID[which(Metadata$Clusters == "malignant")]    ##  Defining Tumor cells subset using Metadata
NormalCells = Metadata$SampleID[which(Metadata$Clusters == "Normal")]       ##  Defining Normal cells subset using Metadata


### Effect size function
### Function Input: data represents the scores for all cells computed by a scoring method, TumorCells vector contains the IDs of all tumor cells, and NormalCells vector contains the IDs of all Normal cells. These IDs would divide the scoring data into tumor and normal cases inside the function

EffectSizeCalculation <- function(data,TumorCells,NormalCells){ 

data_Cancer = data[which(rownames(data) %in% TumorCells),]            ##  Defining Tumor cells scores subset using Tumor cell IDs  
data_Normal = data[which(rownames(data) %in% NormalCells),]           ##  Defining Normal cells scores subset using Normal cell IDs  


OverallES = rbind()

for(j in 1:ncol(data)){
out <- tryCatch(
{
print(j)
effectval = data.frame(colnames(data)[j],cohens_d(data_Cancer[,j],data_Normal[,j]))    ## cohens_d function to compute the effectsize
colnames(effectval)[1] = 'Process'
OverallES = rbind(OverallES,effectval)
},

error = function(cond) {
OverallES = data.frame(Process = "NA", Cohens_d = "NA",CI = "NA","CI_low" = "NA","CI_high" = "NA")

}
)
}
return(OverallES)
}

############################## Effect Size Function Calling for each Signature Scoring Method #########################################

####################### JASMINE  ###############

data = readRDS('JASscores_IDHAstrocytoma.RDS')                     ### Reading combined scoring file of JASMINE for all gene sets
Mean_HNSC_JAS = EffectSizeCalculation(data,TumorCells,NormalCells) ###  Calling effectsize function using JASMINE scores 
colnames(Mean_HNSC_JAS)[2] = 'ES_JAS'
write.table(Mean_HNSC_JAS,'ES_JAS.txt',sep='\t',quote=FALSE,row.names=FALSE)   ### saving the effectsize result for JASMINE scores


####################### ssgsea  ###############

data = readRDS('ssGSEA_IDHAstrocytoma.RDS')								### Reading combined scoring file of ssGSEA for all gene sets
Mean_HNSC_ssgsea = EffectSizeCalculation(data,TumorCells,NormalCells)  ###  Calling effectsize function using ssGSEA scores 
colnames(Mean_HNSC_ssgsea)[2] = 'ES_ssgsea'
write.table(Mean_HNSC_ssgsea,'ES_ssgsea.txt',sep='\t',quote=FALSE,row.names=FALSE)   ### saving the effectsize result for ssGSEA scores

####################### SCSE  ###############

data = readRDS('SCSE_IDHAstrocytoma.RDS')								### Reading combined scoring file of SCSE for all gene sets
Mean_HNSC_SCSE = EffectSizeCalculation(data,TumorCells,NormalCells)			###  Calling effectsize function using SCSE scores 
colnames(Mean_HNSC_SCSE)[2] = 'ES_SCSE'
write.table(Mean_HNSC_SCSE,'ES_SCSE.txt',sep='\t',quote=FALSE,row.names=FALSE)		### saving the effectsize result for SCSE scores

####################### AUCell  ###############

data = readRDS('AUCell_IDHAstrocytoma.RDS')									### Reading combined scoring file of AUCell for all gene sets
Mean_HNSC_AUCell = EffectSizeCalculation(data,TumorCells,NormalCells)		###  Calling effectsize function using AUCell scores 
colnames(Mean_HNSC_AUCell)[2] = 'ES_AUCell'
write.table(Mean_HNSC_AUCell,'ES_AUCell.txt',sep='\t',quote=FALSE,row.names=FALSE)		### saving the effectsize result for AUCell scores


####################### GSVA  ###############

data = readRDS('GSVA_IDHAstrocytoma.RDS')									### Reading combined scoring file of GSVA for all gene sets
Mean_HNSC_GSVA = EffectSizeCalculation(data,TumorCells,NormalCells)			###  Calling effectsize function using GSVA scores 
colnames(Mean_HNSC_GSVA)[2] = 'ES_AUCell'
write.table(Mean_HNSC_GSVA,'ES_GSVA.txt',sep='\t',quote=FALSE,row.names=FALSE)		### saving the effectsize result for GSVA scores


################  Combining the effectsizes for all scoring methods in one data frame #############

Comball = merge(Mean_HNSC_AUCell,Mean_HNSC_JAS,by='Process')
Comball = merge(Comball,Mean_HNSC_SCSE,by='Process')
Comball = merge(Comball,Mean_HNSC_ssgsea,by='Process')
Comball = merge(Comball,Mean_HNSC_GSVA,by='Process')

GSsize = read.table('GeneSetSizes_All_Msigdb_1Feb2021.txt',sep='\t',head=T)   ### Reading the gene set sizes for C2, C3 and Hallmark gene sets
Comball = merge(GSsize,Comball,by='Process')                                  ### Combining the gene set sizes with the ES scores
write.table(Comball,'ES_AllMethods.txt',sep='\t',quote=FALSE,row.names=FALSE)   ### Saving the dataframe containing gene set name, ES scores for 5 scoring methods and Gene set sizes


