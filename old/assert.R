library(phangorn)
options("scipen"=100, "digits"=4)
i=10000000


#read full alignment
data = read.table("46glossesfull.csv", header=TRUE, row.names=1, colClasses='character')
full.alignment = as.matrix(data)
tab <- apply(full.alignment, 2, function(x){table(factor(x, levels=c(0,1)))})
indsToUse <- apply(tab, 2, function(x){x[1] > 1 && x[2] > 1})
full.alignment <- full.alignment[,indsToUse]

#read filtrered alignment
fileName <- paste("execution_results/current_opt_data_large_step", i, ".nex", sep="")
conn <- file(fileName,open="r")
linn <-readLines(conn)
data <- c()
taxanames <- c()
for (j in 7:(length(linn)-3)){
    w = unlist(strsplit(linn[j], split = " "))
    taxanames = c(taxanames, w[5])
    l = unlist(strsplit(linn[j], split = ""))
    data=c(data, l[18:length(l)])
    ncol=length(l)-17
    nrow=length(linn)-9
}
close(conn)
m=matrix(data, nrow, ncol, byrow=TRUE)
rownames(m) <- taxanames

attrNames= read.table(file=paste("execution_results/selected_features_step", i, ".txt", sep=""))
attrNamesVec = c()
for (j in 1:length(attrNames)){
    attrNamesVec = c(attrNamesVec, attrNames[j, 2])
}

colnames(m) <- attrNamesVec
alignment1=as.phyDat(m, type="USER", levels = c("0", "1"))

#read index file
indexes=unlist(read.table(file=paste("execution_results/current_opt_indexes_large_step", i, ".csv", sep=""), sep=","))

#filter full alignment by index file
filtered.data = full.alignment[,indexes]
alignment2 = as.phyDat(filtered.data, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)

#compare with filtered alignment
print(alignment1)
print(alignment2)

#read nj tree
tree1=read.tree(file=paste("execution_results/current_opt_tree_large_step", i, ".nwk", sep=""))

#compute nj tree from filtered alignment
dm <- dist.hamming(alignment1)
write.table(as.matrix(dm), file="dm", sep=" ", col.names=F, row.names=F, quote=F)
tree2 <- NJ(dm)



dm2 <- dist.hamming(alignment2)
tree3 <- NJ(dm2)
#compare trees
d=treedist(tree1, tree2, check.labels=TRUE)
print(paste("***** d: ", d[1] , sep="") )
d2=treedist(tree2, tree3, check.labels=TRUE)
print(paste("***** d2: ", d2[1] , sep="") )
