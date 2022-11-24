library(phangorn)
rm(list=ls())
##set.seed(12383821) ## change this for further analyses

cog.tree=read.tree("IE2011_Cognates_rel_ANNOT.nwk")
data = read.table("46glossesfull.csv", header=TRUE, row.names=1, colClasses='character')
cog.tree$tip.label[cog.tree$tip.label == "Ptg-E"] = "PtgE"
datamat = as.matrix(data)

tab <- apply(datamat, 2, function(x){table(factor(x, levels=c(0,1)))})
indsToUse <- apply(tab, 2, function(x){x[1] > 1 && x[2] > 1})
datamat <- datamat[,indsToUse]
initialdatamat <- datamat
colnames(datamat)
