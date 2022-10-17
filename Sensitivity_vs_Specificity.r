

#############################################################################################################################
######################### SENSITIVITY    VERSUS SPECIFICITY CALCULATION ###############################################
#############################################################################################################################
Dataset = read.table('ConsensusSummaryFile.txt',sep='\t',head=T)

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
