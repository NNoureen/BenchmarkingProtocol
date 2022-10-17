
Data = read.table('CommonCases_ES_IDHAstrocytoma_17March2021.txt',sep='\t',head=T)


#############################################################################################################################
######################### TP, TN, FP and FN calculation for all methods ###############################################
#############################################################################################################################

ConsensusSummary <- function(Data){

######## JASMINE ##################
JAS_cases_common = c('AUCell vs JAS','SCSE vs JAS','AUCell & SCSE & JAS')

JAS_Common_Up = sum(Data$Common[which(Data$Methods %in% JAS_cases_common & Data$ESDir == "Pos")])
JAS_NC_Up = Data$Common[which(Data$Methods == "JAS" & Data$ESDir == "Pos")]

JAS_Up = data.frame("JAS",JAS_Common_Up,JAS_NC_Up,"Pos")
colnames(JAS_Up) = c("Method","Common","NC","ES")

JAS_Common_DN = sum(Data$Common[which(Data$Methods %in% JAS_cases_common & Data$ESDir == "Neg")])
JAS_NC_DN = Data$Common[which(Data$Methods == "JAS" & Data$ESDir == "Neg")]


JAS_Dn = data.frame("JAS",JAS_Common_DN,JAS_NC_DN,"Neg")
colnames(JAS_Dn) = c("Method","Common","NC","ES")


Not_JAS_Common_Up =  Data$Common[which(Data$Methods == "AUCell vs SCSE" & Data$ESDir == "Pos")]
Not_JAS_NC_Up =  sum(Data$Common[which(Data$Methods == c("AUCell","SCSE") & Data$ESDir == "Pos")])


Not_JAS_Up = data.frame("!JAS",Not_JAS_Common_Up,Not_JAS_NC_Up,"Pos")
colnames(Not_JAS_Up) = c("Method","Common","NC","ES")


Not_JAS_Common_Dn =  Data$Common[which(Data$Methods == "AUCell vs SCSE" & Data$ESDir == "Neg")]
Not_JAS_NC_Dn =  sum(Data$Common[which(Data$Methods == c("AUCell","SCSE") & Data$ESDir == "Neg")])

Not_JAS_Dn = data.frame("!JAS",Not_JAS_Common_Dn,Not_JAS_NC_Dn,"Neg")
colnames(Not_JAS_Dn) = c("Method","Common","NC","ES")

######## SCSE ###########################

SCSE_cases_common = c('AUCell vs SCSE','SCSE vs JAS','AUCell & SCSE & JAS')

SCSE_Common_Up = sum(Data$Common[which(Data$Methods %in% SCSE_cases_common & Data$ESDir == "Pos")])
SCSE_NC_Up = Data$Common[which(Data$Methods == "SCSE" & Data$ESDir == "Pos")]

SCSE_Up = data.frame("SCSE",SCSE_Common_Up,SCSE_NC_Up,"Pos")
colnames(SCSE_Up) = c("Method","Common","NC","ES")

SCSE_Common_DN = sum(Data$Common[which(Data$Methods %in% SCSE_cases_common & Data$ESDir == "Neg")])
SCSE_NC_DN = Data$Common[which(Data$Methods == "SCSE" & Data$ESDir == "Neg")]

SCSE_Dn = data.frame("SCSE",SCSE_Common_DN,SCSE_NC_DN,"Neg")
colnames(SCSE_Dn) = c("Method","Common","NC","ES")

Not_SCSE_Common_Up =  Data$Common[which(Data$Methods == "AUCell vs JAS" & Data$ESDir == "Pos")]
Not_SCSE_NC_Up =  sum(Data$Common[which(Data$Methods == c("AUCell","JAS") & Data$ESDir == "Pos")])


Not_SCSE_Up = data.frame("!SCSE",Not_SCSE_Common_Up,Not_SCSE_NC_Up,"Pos")
colnames(Not_SCSE_Up) = c("Method","Common","NC","ES")


Not_SCSE_Common_Dn =  Data$Common[which(Data$Methods == "AUCell vs JAS" & Data$ESDir == "Neg")]
Not_SCSE_NC_Dn =  sum(Data$Common[which(Data$Methods == c("AUCell","JAS") & Data$ESDir == "Neg")])

Not_SCSE_Dn = data.frame("!SCSE",Not_SCSE_Common_Dn,Not_SCSE_NC_Dn,"Neg")
colnames(Not_SCSE_Dn) = c("Method","Common","NC","ES")


######## AUCell ###########################


AUCell_cases_common = c('AUCell vs SCSE','AUCell vs JAS','AUCell & SCSE & JAS')

AUCell_Common_Up = sum(Data$Common[which(Data$Methods %in% AUCell_cases_common & Data$ESDir == "Pos")])
AUCell_NC_Up = Data$Common[which(Data$Methods == "AUCell" & Data$ESDir == "Pos")]

AUCell_Up = data.frame("AUCell",AUCell_Common_Up,AUCell_NC_Up,"Pos")
colnames(AUCell_Up) = c("Method","Common","NC","ES")


AUCell_Common_DN = sum(Data$Common[which(Data$Methods %in% AUCell_cases_common & Data$ESDir == "Neg")])
AUCell_NC_DN = Data$Common[which(Data$Methods == "AUCell" & Data$ESDir == "Neg")]


AUCell_Dn = data.frame("AUCell",AUCell_Common_DN,AUCell_NC_DN,"Neg")
colnames(AUCell_Dn) = c("Method","Common","NC","ES")



Not_AUCell_Common_Up =  Data$Common[which(Data$Methods == "SCSE vs JAS" & Data$ESDir == "Pos")]
Not_AUCell_NC_Up =  sum(Data$Common[which(Data$Methods == c("SCSE","JAS") & Data$ESDir == "Pos")])


Not_AUCell_Up = data.frame("!AUCell",Not_AUCell_Common_Up,Not_AUCell_NC_Up,"Pos")
colnames(Not_AUCell_Up) = c("Method","Common","NC","ES")


Not_AUCell_Common_Dn =  Data$Common[which(Data$Methods == "SCSE vs JAS" & Data$ESDir == "Neg")]
Not_AUCell_NC_Dn =  sum(Data$Common[which(Data$Methods == c("SCSE","JAS") & Data$ESDir == "Neg")])


Not_AUCell_Dn = data.frame("!AUCell",Not_AUCell_Common_Dn,Not_AUCell_NC_Dn,"Neg")
colnames(Not_AUCell_Dn) = c("Method","Common","NC","ES")



CombinedSummary = rbind(JAS_Up,Not_JAS_Up)
CombinedSummary = rbind(CombinedSummary,JAS_Dn)
CombinedSummary = rbind(CombinedSummary,Not_JAS_Dn)
CombinedSummary = rbind(CombinedSummary,SCSE_Up)
CombinedSummary = rbind(CombinedSummary,Not_SCSE_Up)
CombinedSummary = rbind(CombinedSummary,SCSE_Dn)
CombinedSummary = rbind(CombinedSummary,Not_SCSE_Dn)
CombinedSummary = rbind(CombinedSummary,AUCell_Up)
CombinedSummary = rbind(CombinedSummary,Not_AUCell_Up)
CombinedSummary = rbind(CombinedSummary,AUCell_Dn)
CombinedSummary = rbind(CombinedSummary,Not_AUCell_Dn)



return(CombinedSummary)
}

#############################################################################################################################
######################### SENSITIVITY    VERSUS SPECIFICITY CALCULATION ###############################################
#############################################################################################################################
Dataset = ConsensusSummary(Data)

Datacomball = rbind()

for(i in seq(1,nrow(Dataset),2)){
print(i)
subset = Dataset[c(i,i+1),]
sensi = subset$Common[1]/(subset$Common[1]+subset$Common[2])
speci = subset$NC[2]/(subset$NC[1]+subset$NC[2])
accuracy = (subset$Common[1] + subset$NC[2])/ (subset$Common[1]+subset$Common[2]+subset$NC[1]+subset$NC[2])

subsetRes = data.frame(subset$Method[1],subset$ES[1],sensi,speci,accuracy)
colnames(subsetRes) = c('Method','ES','Sensitivity','Specificity','Accuracy')

Datacomball = rbind(Datacomball,subsetRes)

}

write.table(Datacomball,'Specificity_Sensitivity.txt',sep='\t',quote=FALSE,row.names=FALSE)
