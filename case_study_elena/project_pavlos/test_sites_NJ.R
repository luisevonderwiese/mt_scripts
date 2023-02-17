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

colnames(initialdatamat)

phy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)

str(phy)

dm  <- dist.hamming(phy)
treeNJ <- NJ(dm)

#plot(treeNJ)


init.n <- 200 ## the number of sites

max.tries <- 300000
max.fails <- max.tries
old.keep <- 0.9

optimizeSetNJ <- function(data, cog.tree=cog.tree, init.n=ncol(data), max.tries=100000, max.fails=max.tries, minimum.sites=20, maximum.sites=400, old.keep=0.9, proposals=NULL){
    ## iterate to get a distribution to have an idea about the
    ## RF distance
    fails <- 0
    #Matrix max.tries x 2
    proposals <- matrix(NA, ncol=2, nrow=max.tries)
    ## to get an idea of the random distances
    rand.n = 500
    randomRF = vector("numeric", length=rand.n)
    i <- 1
    for(i in 1:rand.n){
        # vector mit n zufälligen zahlen zwischen 1 und #spalten in initialdatamat#
        random.col <- sample(1:ncol(initialdatamat), init.n, replace=TRUE)
        # aus initial data mat die spalten der zufällig gewählten indices auswählen
        datamat = as.matrix(initialdatamat[,random.col])
        # daraus ein sample erzeugen (enthält subset an sites)
        samplephy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
        
        ##optim.pml(treeNJ, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=T, optRate=T)
        #Hamming-distanz-matrix für die taxa / sprachen
        dm  <- dist.hamming(samplephy)
        #neighbor joining baum
        tr <- NJ(dm)


        ##fit = pml(treeNJ, data=samplephy)
        ##tr = optim.pml(fit, model="Mk", rearrangement="stochastic")
        #abstand zum golden standard
        d=treedist(tr, cog.tree, check.labels=TRUE)
        #merken
        randomRF[i] <- d[1]
    }
    pdf("random_histogram_RF.pdf")
    hist(randomRF) ## gyrw sto 65 mexri peripou to 50
    dev.off()
    #bis hier scheint eine art vorschau zu sein
    cur.indexes = sample(1:ncol(initialdatamat), init.n, replace=FALSE)
    datamat = as.matrix(initialdatamat[,cur.indexes])
    new.indexes = cur.indexes
    i <- 1
    curRf = 1000000
    curRf
    n.col.to.replace = round(init.n/15, 0)
    n.col.to.replace
    i <- 1
    for(i in 1:max.tries){
        #alle 100 wird was rausgeschrieben
        if( (i==1) || (i %% 10000 == 0)){
            attrNames <- sort(gsub("^X", "", colnames(initialdatamat)[new.indexes]))
            write.table(attrNames, file=paste("selected_features_step", i, ".txt", sep=""), quote=F, col.names=F, row.names=F)
            write.table(curRf, file=paste("current_accuracy_step", i, ".txt", sep=""), quote=F, col.names=F, row.names=F)
            samplephy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
            dm <- dist.hamming(samplephy)
            tr <- NJ(dm)

            d <- treedist(tr, cog.tree, check.labels =TRUE)
            write.phyDat(x=samplephy, file=paste("current_opt_data_large_step", i, ".nex", sep=""), format="nexus")
            write.table(new.indexes, file=paste("current_opt_indexes_large_step", i, ".csv", sep=""), sep=",", col.names=F, row.names=F, quote=F)
            write.tree(tr, file=paste("current_opt_tree_large_step", i, ".nwk", sep=""))
        }
        print(i)
        newdatamat = datamat
        ## how many columns to replace
        ##n.col.to.replace = sample((n.col.to.replace-1):(n.col.to.replace+1), 1)##rbinom(n=1, size=ncol(data), prob=n.col.to.replace/ncol(data))
        ##if(n.col.to.replace == 0){ n.col.to.replace=1 }
        
        indexes.to.sample.from <- as.vector(1:ncol(initialdatamat))[-new.indexes]
        ## perc of keeping from the old indexes lets.say oldkeep. This is from a bernoulli process
        # wks old.keep entspricht anteil der zuletzt neu gewählten, aber mindestens so hoch wie zuvor
        old.keep <- max(old.keep, length(new.indexes)/ncol(initialdatamat))
        #irgendwie werden die neuen indices recht kompliziert ausgewählt
        #Ein Teil der Indizes vom letzten besten Versuch werden behalten, ein anderer Teil wird neu gezogen
        repeat{
            cur.indexes.to.keep <- sample(c(FALSE, TRUE), size=length(new.indexes), replace=TRUE, prob=c(1-old.keep, old.keep))
            #Man braucht dadaruch, dass man im vorherigen Schritt indizes verwirft, erwartet length(new.indexes)*(1-old.keep) neue
            #Man sampelt aus indexes.to.sample.from, so ergibt suich wahrscheinlichkeit
            prob1 <- min(length(new.indexes)/length(indexes.to.sample.from)*(1-old.keep), 1)
            new.indexes.to.insert <- sample(c(FALSE,TRUE), size=length(indexes.to.sample.from), prob=c(1-prob1, prob1), replace=TRUE)
            new.indexes.to.propose <- c(new.indexes[cur.indexes.to.keep], indexes.to.sample.from[new.indexes.to.insert])
            if( length(new.indexes.to.propose) > minimum.sites & length(new.indexes.to.propose) < maximum.sites){break}
        }
        newdatamat = as.matrix(initialdatamat)[,new.indexes.to.propose]
        samplephy = as.phyDat(newdatamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)

        dm <- dist.hamming(samplephy)
        tr <- NJ(dm)


        ##optim.pml(treeNJ, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=T, optRate=T)
        ##fit = pml(treeNJ, data=samplephy)
        ##tr = optim.pml(fit, model="Mk", rearrangement="stochastic", control = pml.control(trace = 0))
        d=treedist(tr, cog.tree, check.labels=TRUE)
        if(d[1] <= curRf){
            ##print(paste("found a better matrix... current best: ", curRf, " current size: ", length(new.indexes.to.propose)))
            datamat = newdatamat
            #new indices wird nur bei verbesserung gemerkt
            new.indexes = new.indexes.to.propose
            curRf = d[1]
            fails = 0
        }else{
            fails  <- fails+1
            print(c(fails, d[1]))
            #n.col.to.replace = ceiling( round(init.n/10, 0) * (max.fails - fails)/(max.fails+1))
            #scheint abufangen, wenn man in zu tiefem tal landet
            if(fails > max.fails){
                break
            }
        }
        print(paste("***** d: ", d[1], " Current RF: ", curRf, " current size: ", length(new.indexes), " current proposal: ", length(new.indexes.to.propose), "sizes: ", length(new.indexes), ",", sum(cur.indexes.to.keep), ",", sum(new.indexes.to.insert) , sep="") )
        proposals[i,1] <- curRf
        proposals[i,2] <- d[1]
    }
    return(list(datamatrix=datamat, indexes=new.indexes, rf=curRf, size=length(new.indexes)))
}

results=optimizeSetNJ(data=initialdatamat, cog.tree=cog.tree, init.n=150, max.tries=10000000, old.keep=0.99, proposals=proposals)



png("RFpropRF.png")
colnames(results$proposals) <- c("rf", "proposed.rf")
mylm <- lm(proposed.rf~rf, data=as.data.frame(results$proposals))
plot(results$proposals)
abline(a=mylm$coefficients[1], b=mylm$coefficients[2])
dev.off()

## get names
attrNames <- gsub("^X", "", colnames(data)[results$indexes])
write.table(attrNames, file="selected_features.txt", quote=F, col.names=F, row.names=F)

samplephy = as.phyDat(results$datamatrix, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)

dm <- dist.hamming(samplephy)
tr <- NJ(dm)


d <- treedist(tr, cog.tree, check.labels=TRUE)
write.phyDat(x=samplephy, file="opt_data_large.nex", format="nexus")
write.table(results$indexes, file="opt_indexes_large.csv", sep=",", col.names=F, row.names=F, quote=F)
write.tree(tr, file="opt_tree_large.nwk")
d

#plot(tr.node, show.node.label=TRUE)
#write.table(as_tibble(tr.node))
#tr.new = root(tr.node, node=55)

#write.table(as_tibble(tr.node), file="")

#pdf("trees.pdf", width=20, height=10)
#layout(matrix(1:2, nrow=1))
#plot(cog.tree)
#plot(tr.new)
#dev.off()





## likelihood = vector("numeric", length=ncol(data))
## nr=attr(phy, "nr")
## freqs = c(sum(data == 1)/(sum(data == 1 | data == 0)), sum(data == 0)/(sum(data == 1 | data == 0)))

## getp = function(pattern, n=200, freqs=freqs){
##     xx = pattern
##     likes=pml(cog.tree, xx, model="Mk", bf=freqs)
##     lik = likes$logLik
##     rand.lik <- vector("numeric", n)
##     j <- 1
##     for(j in 1:n){
##         rnd = matrix(sample(as.character(pattern), replace=FALSE), ncol=1)
##         row.names(rnd) = rownames(data)
##         xx = as.phyDat(rnd, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data))
##         likes=pml(cog.tree, xx, model="Mk", bf=freqs)
##         rlik = likes$logLik
##         paste(as.character(xx), as.character(pattern))
##         rand.lik[j] <- rlik
##     }
##     pvalue = (sum( rand.lik >= lik))/n + 1e-6
##     return(pvalue)
## }

## pvals = vector("numeric", nr)

## i <- 1
## as.character(phy[,i])

## for(i in 1:nr){
##     print(i)
##     pvals[i] <- getp(phy[,i], n=1000, freqs=freqs)
## }



## pval.thr=0.4
## good.index = vector(mode="logical",length=ncol(data))
## for(i in 1:ncol(data)){
##     pat.index = attr(phy, which="index")[i]
##     xx = phy[,pat.index]
##     likes=pml(cog.tree, xx, model="Mk", bf=freqs)
##     likelihood[i] = likes$logLik
##     if(pvals[pat.index] < pval.thr){
##         good.index[i] = TRUE
##     }else{
##         good.index[i] = FALSE
##     }
## }

##                                         #plot(density(likelihood))
## #Sys.setenv(DISPLAY="localhost:12.0")
## #sort(likelihood, decreasing=TRUE)
## #order(likelihood, decreasing=TRUE)[1:10]

## #data[,order(likelihood, decreasing=TRUE)[1:10]]
## pvals                                        #likelihood[2] 
## sum(good.index)


## write.table(data[,good.index], file="high_likelihood_exp2_pval0.2.tsv", col.names=FALSE, row.names=TRUE, sep="\t", quote=FALSE)

## #as.phyDat(rnd, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data))
## goodphy = as.phyDat(as.matrix(data[,good.index]), type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data))
## #
## #optim.pml(treeNJ, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=T, optRate=T)
## fit = pml(treeNJ, data=goodphy)
## tr = optim.pml(fit, model="Mk", optGamma=TRUE, rearrangement="stochastic", optEdge=TRUE, optNni=TRUE)




## library(ape)
## library(tidytree)

## bestp = 0

## res <- list()
## for(j in 1:100){
##     print(j)
##     pval.thr=j/100
##     good.index = vector(mode="logical",length=ncol(data))
##     for(i in 1:ncol(data)){
##         pat.index = attr(phy, which="index")[i]
##         ##xx = phy[,pat.index]
##         ##likes=pml(cog.tree, xx, model="Mk", bf=freqs)
##         ##likelihood[i] = likes$logLik
##         if(pvals[pat.index] < pval.thr){
##             good.index[i] = TRUE
##         }else{
##             good.index[i] = FALSE
##         }
##     }
##     goodphy = as.phyDat(as.matrix(data[,good.index]), type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data))
##     ##
##     ##optim.pml(treeNJ, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=T, optRate=T)
##     fit = pml(treeNJ, data=goodphy)
##     tr = optim.pml(fit, model="Mk", optGamma=TRUE, rearrangement="stochastic", optEdge=TRUE, optNni=TRUE)
##     d=treedist(tr$tree, cog.tree, check.labels=TRUE)
##     res[[j]] <- c(pval.thr, d)
## }
## res.mat=matrix(unlist(res), ncol=5, byrow=TRUE)
## plot(res.mat[,1], res.mat[,2])

## #################   get good trees by comparing a random dataset to the target tree #####################


## init.n=200
## optimizeSet <- function(target, data, init.n=100){
##     ## iterate to get a distribution to have an idea about the
##     ## RF distance
##     max.fails = 500
##     fails <- 0
##     rand.n = 100
##     randomRF = vector("numeric", length=rand.n)
##     i <- 1
##     for(i in 1:rand.n){
##         random.col <- sample(1:ncol(data), init.n, replace=TRUE)
##         datamat = as.matrix(data[,random.col])
##         samplephy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
##         ##optim.pml(treeNJ, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=T, optRate=T)
##         fit = pml(treeNJ, data=samplephy)
##         tr = optim.pml(fit, model="Mk", rearrangement="stochastic")
##         d=treedist(tr$tree, cog.tree, check.labels=TRUE)
##         randomRF[i] <- d[1]
##     }

##     cur.indexes = random.col
##     new.indexes = cur.indexes
##     i <- 1
##     curRf = d[1]##randomRF[length(randomRF)] ## the last rf
##     curRf
##     n.col.to.replace = round(init.n/10, 0)
##     n.col.to.replace
##     for(i in 1:1000){
##         print(i)
##         newdatamat = datamat
##         ## how many columns to replace
##         col.to.replace = sample(1:ncol(newdatamat), n.col.to.replace, replace=FALSE)
##         new.col = sample(1:ncol(data), n.col.to.replace, replace=TRUE)
##         newdatamat[,col.to.replace] = as.matrix(data)[,new.col]
##         samplephy = as.phyDat(newdatamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
##         ##optim.pml(treeNJ, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=T, optRate=T)
##         fit = pml(treeNJ, data=samplephy)
##         tr = optim.pml(fit, model="Mk", rearrangement="stochastic", control = pml.control(trace = 0))
##         d=treedist(tr$tree, cog.tree, check.labels=TRUE)
##         if(d[1] <= curRf){
##             print(paste("found a better matrix... current best: ", curRf))
##             datamat = newdatamat
##             cur.indexes[col.to.replace] = new.col
##             new.indexes = cur.indexes
##             curRf = d[1]
##             fails = 0
##         }else{
##             fails  <- fails+1
##             print(c(fails, d[1]))
##             n.col.to.replace = ceiling( round(init.n/10, 0) * (max.fails - fails)/(max.fails+1))
##             if(fails > max.fails){
##                 break
##             }
##         }
##         print(paste("******************************** Current RF: ", curRf, " ncolreplace:", n.col.to.replace))
##     }
## }
