library(MAST)
library(Seurat)
library('ggplot2')
library('gplots')
library('dplyr')
stringsAsFactors=FALSE
library(circlize)
library(stringr)
library(data.table)
library('GSA')
###############################

data = readRDS('IDHAstrocytoma_GE_20210311.RDS')		### Read the Gene expression data 
genes = rownames(data)									### create the list of genes present in the data

##############################

DnGenesList <-  read.table('IDHastrocytoma_MAST_DEGS_15March2021.txt',sep='\t',head=T)		### read the DEGS file generated using MAST function
DnGenesList = DnGenesList[which(DnGenesList$avg_logFC < 0),]								### Shortlist the down regulated genes based on logFC values < 0
NoiseGenesList <- genes[-which(genes %in% DnGenesList)] 									### Take remaining genes as noise list	


GSsizesList = c(50,100,150,200,300)														### Intialize the gene set sizes for gold standard gene sets


###### section to generate gold standard gene sets for each gene set size intialized above

for(j in 1:length(GSsizesList)){
L1 = GSsizesList[j]
print(L1)


###### for each gene set size 200 gene sets would be generated with randomly sampled genes and 4 noise levels

for(k in 1:200){
origid = paste("GS",k,sep="")
origid = paste(origid, L1,sep="_")
originalid = paste("OS",origid,sep="_")

OriginalSet = sample(DnGenesList,L1)
OS = c(originalid,"Random",OriginalSet)


sink("IDHAstrocytoma_DN.gmt",append=T)			### save the DN regulated gold standard data set
writeLines(OS,sep='\t')
writeLines("",sep='\n')
sink()

Noise1 = L1 * 0.2					### 20% Noise Level generation 
Noise2 = L1 * 0.4					### 40% Noise Level generation 
Noise3 = L1 * 0.6					### 60% Noise Level generation 
Noise4 = L1 * 0.8					### 80% Noise Level generation 



NoiseLevels = c(Noise1,Noise2,Noise3,Noise4)

for(i in 1:length(NoiseLevels)){
subsetsize = NoiseLevels[i]
OS_subset = L1 - subsetsize
Noisesubset = c(sample(OriginalSet,OS_subset),sample(NoiseGenesList,subsetsize))
id = paste("N",i,sep="")
id = paste(id,origid,sep="_")
Noiseadd = c(id,"Random",Noisesubset)
sink("IDHAstrocytoma_DN.gmt",append = T)
writeLines(Noiseadd,sep='\t')
writeLines("",sep='\n')
sink()

}
}
}
