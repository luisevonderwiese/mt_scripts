library(phangorn)

data = read.table("46glossesfull.csv", header=TRUE, row.names=1, colClasses='character')
datamat = as.matrix(data)
phy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
dim(datamat)
common <- vector(mode="numeric", length=1000)
for(i in 1:1000){
    s1 <- sample(1:ncol(datamat), 200, replace=F)
    s2 <- sample(1:ncol(datamat), 200, replace=F)
    x1 <- attr(x=phy, "index")[s1]
    x2 <- attr(x=phy, "index")[s2]
    common[i] <- sum(x1 %in% x2)
    print(paste(i, " ", sum(x1 %in% x2), sep=""))
}


for(i in 1:11){
    j <- as.integer(i*10000)
    f <- paste("selected_features_step", j, ".txt", sep="")
    a <- read.table(f)
                                        #a <- read.table("selected_features_step110000.txt")
    b <- read.table("../46_gloses_200sites_uniqSites/selected_features.txt")
    ind1 <- sapply(a[,1], function(x){ which( x == gsub("^X", "", colnames(data) ) )})
    ind2 <- sapply(b[,1], function(x){ which( x == gsub("^X", "", colnames(data) ) )})

    x1 <- attr(x=phy, "index")[ind1]
    x2 <- attr(x=phy, "index")[ind2]
    
    print(paste(j, " ", sum(x1 %in% x2), sep=""))
}



