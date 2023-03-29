# Before starting:
#1. make sure you have installed the "ape" package.
#2. the file code.R is in the working directory.

library("ape")
library(R.utils)
library(tools)
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
        deltaA <- withTimeout({delta(trait,tree,lambda0,se,sim,thin,burn)}, timeout = 10, onTimeout = "warning")
        if (is.null(deltaA)){
          deltaA<--1
          print(trait)
        }
        #deltaA <- delta(trait,tree,lambda0,se,sim,thin,burn)
        print(deltaA)
        deltas = append(deltas, deltaA)
    }
    tree.name.only = file_path_sans_ext(basename(tree.name))
    file.name = paste("deltas/", tree.name.only, "_morpho_filtered_indsToUse.deltas", sep = "")
    write.table(deltas, file=file.name, sep=",", col.names=F, row.names=F, quote=F)


}


#files <- list.files(path="rerooted_trees/geo_duration/inner_nodes/", pattern="*.tree", full.names=TRUE, recursive=FALSE)
#lapply(files, apply.delta.statistics)
apply.delta.statistics("rerooted_trees/cognate_ie_compatible/cognate_ie_compatible_rd.tree")
apply.delta.statistics("rerooted_trees/geo_duration/geo_duration_rd.tree")



