
library('ggplot2')
library('gplots')
library('dplyr')
library('ggpubr')
stringsAsFactors=FALSE
library(stringr)
library(data.table)

###################### Consensus Function Calling ##############

Data = read.table('ES_AllMethods.txt',sep='\t',head=T)         ### Reading ES scores file containing ES for all scoring methods with gene set sizes
Data = Data[which(Data$GSsize > 20),]						### Removing the gene sets with size less than 20
Res = Sensitivityfunc(Data)								### Calling consensus function
Res$Dataset = "IDHAstrocytoma"
write.table(Res,'CommonCases_ES_IDHAstrocytoma_17March2021.txt',sep='\t',quote=FALSE,row.names=FALSE)   ### saving the results

######################  Consensus Function Definition ###########
## Function Input:  ES scores for all scoring method for a single cell data set
## Output:  Combination of Up (Pos) and DN (Neg) cases common across 2 or more methods and for each scoring method . These are used to define TP, TN, FP and FN .

Sensitivityfunc <- function(Data){
CommonCases = rbind()

Case1 = Data[which((Data$ES_AUCell >= 1 & Data$ES_SCSE >= 1 & Data$ES_JAS < 1)),]

subsetcommon = data.frame(nrow(Case1),'AUCell vs SCSE','Pos')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)

Case2 = Data[which((Data$ES_AUCell >= 1 & Data$ES_JAS >= 1 & Data$ES_SCSE < 1)),]

subsetcommon = data.frame(nrow(Case2),'AUCell vs JAS','Pos')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)


Case3 = Data[which((Data$ES_SCSE >= 1 & Data$ES_JAS >= 1 & Data$ES_AUCell < 1)),]

subsetcommon = data.frame(nrow(Case3),'SCSE vs JAS','Pos')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)



Case4 = Data[which((Data$ES_AUCell <= -1 & Data$ES_SCSE <= -1 & Data$ES_JAS > -1)),]

subsetcommon = data.frame(nrow(Case4),'AUCell vs SCSE','Neg')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)



Case5 = Data[which((Data$ES_AUCell <= -1 & Data$ES_JAS <= -1 & Data$ES_SCSE > -1)),]

subsetcommon = data.frame(nrow(Case5),'AUCell vs JAS','Neg')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)

Case6 = Data[which((Data$ES_SCSE <= -1 & Data$ES_JAS <= -1 & Data$ES_AUCell > -1)),]

subsetcommon = data.frame(nrow(Case6),'SCSE vs JAS','Neg')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)


Case7 = Data[which(Data$ES_AUCell >= 1 & Data$ES_SCSE < 1 & Data$ES_JAS < 1),]
subsetcommon = data.frame(nrow(Case7),'AUCell','Pos')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)

Case8 = Data[which(Data$ES_SCSE >= 1 & Data$ES_JAS < 1 & Data$ES_AUCell < 1),]
subsetcommon = data.frame(nrow(Case8),'SCSE','Pos')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)



Case9 = Data[which(Data$ES_JAS >= 1 & Data$ES_SCSE < 1 & Data$ES_AUCell < 1),]
subsetcommon = data.frame(nrow(Case9),'JAS','Pos')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)



Case10 = Data[which(Data$ES_AUCell <= -1 & Data$ES_SCSE > -1 & Data$ES_JAS > -1),]
subsetcommon = data.frame(nrow(Case10),'AUCell','Neg')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)


Case11 = Data[which(Data$ES_SCSE <= -1 & Data$ES_JAS > -1 & Data$ES_AUCell > -1),]
subsetcommon = data.frame(nrow(Case11),'SCSE','Neg')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)


Case12 = Data[which(Data$ES_JAS <= -1 & Data$ES_SCSE > -1 & Data$ES_AUCell > -1),]
subsetcommon = data.frame(nrow(Case12),'JAS','Neg')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)


Case13 = Data[which(Data$ES_JAS >= 1 & Data$ES_SCSE >= 1 & Data$ES_AUCell >= 1),]
subsetcommon = data.frame(nrow(Case13),'AUCell & SCSE & AUCell','Pos')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)

Case14 = Data[which(Data$ES_JAS <= -1 & Data$ES_SCSE <= -1 & Data$ES_AUCell <= -1),]
subsetcommon = data.frame(nrow(Case14),'AUCell & SCSE & AUCell','Neg')
colnames(subsetcommon) = c('Common','Methods','ESDir')
CommonCases = rbind(CommonCases,subsetcommon)

return(CommonCases)
}


