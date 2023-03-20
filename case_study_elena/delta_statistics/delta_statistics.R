# Before starting:
#1. make sure you have installed the "ape" package.
#2. the file code.R is in the working directory.

library("ape")
library(R.utils)
source("code.R")

#SOME PARAMETERS... 
lambda0 <- 0.1   #rate parameter of the proposal 
se      <- 0.5   #standard deviation of the proposal
sim     <- 10000 #number of iterations
thin    <- 10    #we kept only each 10th iterate 
burn    <- 100   #100 iterates are burned-in

apply.delta.statistics <- function(tree.name) {
    data = read.table("46glossesfull.csv", header=TRUE, row.names=1, colClasses='character')
    taxa_names = row.names(data)
    datamat = as.matrix(data)
    tab <- apply(datamat, 2, function(x){table(factor(x, levels=c(0,1)))})
    indsToUse <- apply(tab, 2, function(x){x[1] > 1 && x[2] > 1})
    datamat <- datamat[,indsToUse]
    tree<-read.tree(file=tree.name)
    tree$tip.label[tree$tip.label == "Ptg-E"] = "PtgE"
    deltas = c()
    for (col in 1:ncol(datamat)) {
        trait<-datamat[, col]
        deltaA <- withTimeout({delta(trait,tree,lambda0,se,sim,thin,burn)}, timeout = 30, onTimeout = "warning")
        if (is.null(deltaA)){
          deltaA<--1
        }
        #deltaA <- delta(trait,tree,lambda0,se,sim,thin,burn)
        print(deltaA)
        deltas = append(deltas, deltaA)
    }
    write.table(deltas, file="deltas.csv", sep=",", col.names=F, row.names=F, quote=F)


}

apply.delta.statistics("cognate_ie_compatible.tree.rooted.tree")

