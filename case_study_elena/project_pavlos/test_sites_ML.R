library(phangorn)
rm(list=ls())
##set.seed(12383821) ## change this for further analyses



optimizeSetNJ <- function(data, cog.tree=cog.tree, init.n=ncol(data), max.tries=100000, max.fails=max.tries, minimum.sites=20, maximum.sites=400, old.keep=0.9, start.idx=1, cur.indexes=NULL){
    ## iterate to get a distribution to have an idea about the
    ## RF distance
    fails <- 0
    if(length(cur.indexes)==0){
        cur.indexes = sample(1:ncol(initialdatamat), init.n, replace=FALSE)
    }
    new.indexes = cur.indexes
    datamat = as.matrix(initialdatamat[,cur.indexes])
    samplephy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)

    system("mkdir raxml")
    write.phyDat(x=samplephy, file="raxml/current_alignment.phy", format="phylip")
    system("./raxml-ng --search1 --msa raxml/current_alignment.phy --seed 2 --model BIN+G --blopt nr_safe --prefix raxml/")
    bestTree=read.tree("raxml/.raxml.bestTree")
    system("rm -rf raxml")

    d=treedist(bestTree, cog.tree, check.labels=TRUE)
    curRf = d[1]

    n.col.to.replace = round(init.n/15, 0)
    n.col.to.replace

    dir.create(file.path("ml_results/"), showWarnings = FALSE)
    for(i in start.idx:max.tries){
        #alle 100 wird was rausgeschrieben
        if( (i==1) || (i %% 10000 == 0)){
            attrNames <- sort(gsub("^X", "", colnames(initialdatamat)[new.indexes]))
            samplephy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
            write.table(attrNames, file=paste("ml_results/ml_selected_features_step", i, ".txt", sep=""), quote=F, col.names=F, row.names=F)
            write.table(curRf, file=paste("ml_results/ml_current_accuracy_step", i, ".txt", sep=""), quote=F, col.names=F, row.names=F)
            write.phyDat(x=samplephy, file=paste("ml_results/ml_current_opt_data_large_step", i, ".nex", sep=""), format="nexus")
            write.table(new.indexes, file=paste("ml_results/ml_current_opt_indexes_large_step", i, ".csv", sep=""), sep=",", col.names=F, row.names=F, quote=F)
            write.tree(bestTree, file=paste("ml_results/ml_current_opt_tree_large_step", i, ".nwk", sep=""))
        }
        print(i)
        newdatamat = datamat
        indexes.to.sample.from <- as.vector(1:ncol(initialdatamat))[-new.indexes]
        old.keep <- max(old.keep, length(new.indexes)/ncol(initialdatamat))
        repeat{
            cur.indexes.to.keep <- sample(c(FALSE, TRUE), size=length(new.indexes), replace=TRUE, prob=c(1-old.keep, old.keep))
            prob1 <- min(length(new.indexes)/length(indexes.to.sample.from)*(1-old.keep), 1)
            new.indexes.to.insert <- sample(c(FALSE,TRUE), size=length(indexes.to.sample.from), prob=c(1-prob1, prob1), replace=TRUE)
            new.indexes.to.propose <- c(new.indexes[cur.indexes.to.keep], indexes.to.sample.from[new.indexes.to.insert])
            if( length(new.indexes.to.propose) > minimum.sites & length(new.indexes.to.propose) < maximum.sites){break}
        }
        newdatamat = as.matrix(initialdatamat)[,new.indexes.to.propose]
        samplephy = as.phyDat(newdatamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)

        system("mkdir raxml")
        write.phyDat(x=samplephy, file="raxml/current_alignment.phy", format="phylip")
        system("./raxml-ng --search1 --msa raxml/current_alignment.phy --seed 2 --model BIN+G --blopt nr_safe --prefix raxml/")
        tr=read.tree("raxml/.raxml.bestTree")
        system("rm -rf raxml")

        d=treedist(tr, cog.tree, check.labels=TRUE)
        if(d[1] <= curRf){
            datamat = newdatamat
            new.indexes = new.indexes.to.propose
            curRf = d[1]
            bestTree = tr
            fails = 0
        }else{
            fails  <- fails+1
            print(c(fails, d[1]))
            if(fails > max.fails){
                break
            }
        }
        print(paste("***** d: ", d[1], " Current RF: ", curRf, " current size: ", length(new.indexes), " current proposal: ", length(new.indexes.to.propose), "sizes: ", length(new.indexes), ",", sum(cur.indexes.to.keep), ",", sum(new.indexes.to.insert) , sep="") )
    }
    return(list(datamatrix=datamat, indexes=new.indexes, rf=curRf, size=length(new.indexes)))
}

cog.tree=read.tree("IE2011_Cognates_rel_ANNOT.nwk")
data = read.table("46glossesfull.csv", header=TRUE, row.names=1, colClasses='character')
cog.tree$tip.label[cog.tree$tip.label == "Ptg-E"] = "PtgE"
datamat = as.matrix(data)

tab <- apply(datamat, 2, function(x){table(factor(x, levels=c(0,1)))})
indsToUse <- apply(tab, 2, function(x){x[1] > 1 && x[2] > 1})
datamat <- datamat[,indsToUse]
initialdatamat <- datamat
colnames(datamat)
colnames(initialdatamat)

i = 1
indexes = NULL
#indexes = unlist(read.table(file=paste("ml_results/ml_current_opt_indexes_large_step", i, ".csv", sep=""), sep=","))
print(indexes)
results=optimizeSetNJ(data=initialdatamat, cog.tree=cog.tree, init.n=150, max.tries=10000000, old.keep=0.99, start.idx=i, cur.indexes=indexes)


