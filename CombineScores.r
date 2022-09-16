library('dplyr')
stringsAsFactors=FALSE
library(stringr)
library(data.table)
library(effectsize)

######################### Combining JASMINE scores for C2, C3 and Hallmark gene sets in one object ###########################


JAS1 =  as.data.frame(fread("JAS_IDHA_HallmarksResult.txt"))            ### Reading Hallmark gene sets scores computed by JASMINE function
JAS2 = as.data.frame(fread("JAS_IDHA_C2_GSets.txt"))                    ### Reading C2 gene sets scores computed by JASMINE function
JAS3 = as.data.frame(fread("JAS_IDHA_C3_GSets.txt"))                    ### Reading C3 gene sets scores computed by JASMINE function

HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)                           ### Combining the scores 
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb                        
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'JASscores_IDHAstrocytoma.RDS')                 ### Saving the combined JASMINE scores in an R object 



######################### Combining AUCell scores for C2, C3 and Hallmark gene sets in one object ##################################################

JAS1 =  as.data.frame(fread("AUCell_IDHA_HallmarksResult.txt"))       ### Reading Hallmark gene sets scores computed by AUCell function
JAS2 = as.data.frame(fread("AUCell_IDHA_C2.txt"))                     ### Reading C2 gene sets scores computed by AUCell function
JAS3 = as.data.frame(fread("AUCell_IDHA_C3_GSets.txt"))				  ### Reading C3 gene sets scores computed by AUCell function


HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)                          ### Combining the scores 
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'AUCell_IDHAstrocytoma.RDS')                 ### Saving the combined AUCell scores in an R object 



######################### Combining SCSE scores for C2, C3 and Hallmark gene sets in one object##################################################

JAS1 =  as.data.frame(fread("SCSE_IDHA_HallmarksResult.txt"))               ### Reading Hallmark gene sets scores computed by SCSE function
JAS2 = as.data.frame(fread("SCSE_IDHA_C2_GSets.txt"))                       ### Reading C2 gene sets scores computed by SCSE function
JAS3 = as.data.frame(fread("SCSE_IDHA_C3.txt"))                             ### Reading C3 gene sets scores computed by SCSE function


HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)                               ### Combining the scores 
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'SCSE_IDHAstrocytoma.RDS')                    ### Saving the combined SCSE scores in an R object 

########################## Combining GSVA scores for C2, C3 and Hallmark gene sets in one object#################################################

JAS1 =  as.data.frame(fread("GSVA_IDHA_HallmarksResult.txt"))                ### Reading Hallmark gene sets scores computed by GSVA function
JAS2 = as.data.frame(fread("GSVA_IDHA_C2_GSets.txt"))						 ### Reading C2 gene sets scores computed by GSVA function
JAS3 = as.data.frame(fread("GSVA_IDHA_C3.txt"))                              ### Reading C3 gene sets scores computed by GSVA function


HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)                                ### Combining the scores 
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'GSVA_IDHAstrocytoma.RDS')                    ### Saving the combined GSVA scores in an R object 

##########################Combining ssGSEA scores for C2, C3 and Hallmark gene sets in one object#################################################

JAS1 =  as.data.frame(fread("ssGSEA_IDHA_HallmarksResult.txt"))                ### Reading Hallmark gene sets scores computed by ssGSEA function
JAS2 = as.data.frame(fread("ssGSEA_IDHA_C2_GSets.txt"))						   ### Reading C2 gene sets scores computed by ssGSEA function
JAS3 = as.data.frame(fread("ssGSEA_IDHA_C3.txt"))							   ### Reading C3 gene sets scores computed by ssGSEA function


HNSC_C5_Results = rbind(JAS1,JAS2)
HNSC_C5_Results = rbind(HNSC_C5_Results,JAS3)                                  ### Combining the scores 
row.names(HNSC_C5_Results) = HNSC_C5_Results$Msigb
lastidx = ncol(HNSC_C5_Results)
HNSC_C5_Results = HNSC_C5_Results[,-lastidx]
HNSC_C5_Results = t(HNSC_C5_Results)
saveRDS(HNSC_C5_Results,'ssGSEA_IDHAstrocytoma.RDS')                  ### Saving the combined ssGSEA scores in an R object 

