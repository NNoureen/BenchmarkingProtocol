library('dplyr')
stringsAsFactors=FALSE
library(stringr)
library(data.table)
library(effectsize)

######################### Combining JASMINE scores for All gene sets ###########################


JAS1 =  as.data.frame(fread("JAS_IDHA_HallmarksResult.txt"))
JAS2 = as.data.frame(fread("JAS_IDHA_C2_GSets.txt"))
JAS3 = as.data.frame(fread("JAS_IDHA_C3_GSets.txt"))

HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'JASscores_IDHAstrocytoma.RDS')



######################### Combining AUCell scores for All gene sets##################################################

JAS1 =  as.data.frame(fread("AUCell_IDHA_HallmarksResult.txt"))
JAS2 = as.data.frame(fread("AUCell_IDHA_C3.txt"))
JAS3 = as.data.frame(fread("AUCell_IDHA_C2_GSets.txt"))


HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'AUCell_IDHAstrocytoma.RDS')



######################### SCSE ##################################################

JAS1 =  as.data.frame(fread("SCSE_IDHA_HallmarksResult.txt"))
JAS2 = as.data.frame(fread("SCSE_IDHA_C2_GSets.txt"))
JAS3 = as.data.frame(fread("SCSE_IDHA_C3.txt"))


HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'SCSE_IDHAstrocytoma.RDS')

########################## GSVA #################################################

JAS1 =  as.data.frame(fread("GSVA_IDHA_HallmarksResult.txt"))
JAS2 = as.data.frame(fread("GSVA_IDHA_C2_GSets.txt"))
JAS3 = as.data.frame(fread("GSVA_IDHA_C3.txt"))


HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'GSVA_IDHAstrocytoma.RDS')

########################## ssGSEA #################################################

JAS1 =  as.data.frame(fread("ssGSEA_IDHA_HallmarksResult.txt"))
JAS2 = as.data.frame(fread("ssGSEA_IDHA_C2_GSets.txt"))
JAS3 = as.data.frame(fread("ssGSEA_IDHA_C3.txt"))


HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'ssGSEA_IDHAstrocytoma.RDS')



############################## Effect Size Function #########################################
HNSC_Classifications = read.table('IDHAstrocytoma_Metadata.txt',sep='\t',head=T)
HNSC_Stem = HNSC_Classifications$SampleID[which(HNSC_Classifications$Clusters == "Tumor")]
HNSC_Diff = HNSC_Classifications$SampleID[which(HNSC_Classifications$Clusters == "Normal")]


EffectSizeCalculation <- function(data,HNSC_Stem,HNSC_Diff){

data_Cancer = data[which(rownames(data) %in% HNSC_Stem),]
data_Normal = data[which(rownames(data) %in% HNSC_Diff),]


OverallES = rbind()

for(j in 1:ncol(data)){
out <- tryCatch(
{
print(j)
effectval = data.frame(colnames(data)[j],cohens_d(data_Cancer[,j],data_Normal[,j]))
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

############################## Effect Size Function Calling for each Method #########################################

####################### JASMINE Cohens_d ###############

data = readRDS('JASscores_IDHAstrocytoma.RDS')
Mean_HNSC_JAS = EffectSizeCalculation(data,HNSC_Stem,HNSC_Diff)
colnames(Mean_HNSC_JAS)[2] = 'ES_JAS'
write.table(Mean_HNSC_JAS,'ES_JAS.txt',sep='\t',quote=FALSE,row.names=FALSE)


####################### ssgsea Cohens_d ###############

data = readRDS('ssGSEA_IDHAstrocytoma.RDS')
Mean_HNSC_ssgsea = EffectSizeCalculation(data,HNSC_Stem,HNSC_Diff)
colnames(Mean_HNSC_ssgsea)[2] = 'ES_ssgsea'
write.table(Mean_HNSC_ssgsea,'ES_ssgsea.txt',sep='\t',quote=FALSE,row.names=FALSE)

####################### SCSE Cohens_d ###############

data = readRDS('SCSE_IDHAstrocytoma.RDS')
Mean_HNSC_SCSE = EffectSizeCalculation(data,HNSC_Stem,HNSC_Diff)
colnames(Mean_HNSC_SCSE)[2] = 'ES_SCSE'
write.table(Mean_HNSC_SCSE,'ES_SCSE.txt',sep='\t',quote=FALSE,row.names=FALSE)

####################### AUCell Cohens_d ###############

data = readRDS('AUCell_IDHAstrocytoma.RDS')
Mean_HNSC_AUCell = EffectSizeCalculation(data,HNSC_Stem,HNSC_Diff)
colnames(Mean_HNSC_AUCell)[2] = 'ES_AUCell'
write.table(Mean_HNSC_AUCell,'ES_AUCell.txt',sep='\t',quote=FALSE,row.names=FALSE)


####################### GSVA Cohens_d ###############

data = readRDS('GSVA_IDHAstrocytoma.RDS')
Mean_HNSC_GSVA = EffectSizeCalculation(data,HNSC_Stem,HNSC_Diff)
colnames(Mean_HNSC_GSVA)[2] = 'ES_AUCell'
write.table(Mean_HNSC_GSVA,'ES_GSVA.txt',sep='\t',quote=FALSE,row.names=FALSE)


################

Comball = merge(Mean_HNSC_AUCell,Mean_HNSC_JAS,by='Process')
Comball = merge(Comball,Mean_HNSC_SCSE,by='Process')
Comball = merge(Comball,Mean_HNSC_ssgsea,by='Process')
Comball = merge(Comball,Mean_HNSC_GSVA,by='Process')

write.table(Comball,'ES_AllMethods.txt',sep='\t',quote=FALSE,row.names=FALSE)


