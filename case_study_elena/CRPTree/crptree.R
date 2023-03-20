library(CRPTree)

update_time<-function(tree, coal_times){
    #coal_times<-cumsum(coalescent.intervals(tree)$interval.length)
    #tiplabels();nodelabels()
    ##sort before update to avoid problems
    #xx1<-sort(n.t,index.return=T)
    old.edge<-tree$edge
    n.sample <- tree$Nnode + 1
    t.tot <- max(ape::node.depth.edgelength(tree))
    n.t <- t.tot - ape::node.depth.edgelength(tree) ##gives the node length
    new.n.t<-n.t
    #order nodes according to length, then coalescent times are in reverse order
    xx1<-sort(new.n.t,index.return=T)
    index<-(2*n.sample-1):(n.sample+1)
    for (j in (n.sample+1):(2*n.sample-1)){
        old.edge[which(tree$edge[,1]==xx1$ix[j]),1]<-index[j-n.sample]
        old.edge[which(tree$edge[,2]==xx1$ix[j]),2]<-index[j-n.sample]
    }
    new.n.t[(n.sample+1):(2*n.sample-1)]<-rev(coal_times)
    #If we sort them, we can get the correspondence between leaves and coal. times
    xx<-sort(new.n.t,index.return=T)
    new.edge<-old.edge
    for (j in (n.sample+1):(2*n.sample-1)){
        new.edge[which(old.edge[,1]==xx$ix[j]),1]<-index[j-n.sample]
        new.edge[which(old.edge[,2]==xx$ix[j]),2]<-index[j-n.sample]
    }
    new.edge.length<-new.n.t[new.edge[,1]]-new.n.t[new.edge[,2]]
    new.tree<-tree
    new.tree$edge<-new.edge
    new.tree$edge.length<-new.edge.length
    #new.tree$tip.label<-rev(tree$tip.label)
    #tree2<-write.tree(new.tree)
    #new.tree2<-read.tree(text=tree2)
    trees <- c(tree, new.tree)
    trees <- .compressTipLabel(trees)
    t2 <- trees[[2]]
    return(t2)
}

rank.tree <- function(tree) {
        class(tree) <- 'phylo'
        inter_coal_times<-coalescent.intervals(tree)$interval.length
        #m = min(inter_coal_times)
        #inter_coal_times<-inter_coal_times-m
        inter_coal_times[which(inter_coal_times<=0.001)]<-.1
        coal_times<-cumsum(inter_coal_times)
        print(coal_times)
        tree2<-update_time(tree,coal_times)
        tree2<-process_tree(tree2)
        return(tree2)
}


apply.crp.method <- function(tree.name) {
    data = read.table("46glossesfull.csv", header=TRUE, row.names=1, colClasses='character')
    taxa_names = row.names(data)
    datamat = as.matrix(data)
    tab <- apply(datamat, 2, function(x){table(factor(x, levels=c(0,1)))})
    indsToUse <- apply(tab, 2, function(x){x[1] > 1 && x[2] > 1})
    datamat <- datamat[,indsToUse]
    tree<-read.tree(file=tree.name)
    tree$tip.label[tree$tip.label == "Ptg-E"] = "PtgE"
    tree = rank.tree(tree)
    for (col in 1:ncol(datamat)) {
        trait<-datamat[, col]
        converted.trait <- c()
        for (value in trait) {
            if (value == "0"){
                converted.trait = append(converted.trait, -1)
            } else if (value == "1"){
                converted.trait = append(converted.trait, -2)
            } else {
                converted.trait = append(converted.trait, -1)
            }
        }
        cur.tree<-list(edge = tree$edge, tip.label = converted.trait, Nnode = tree$Nnode)
        cur.tree<-rank.tree(cur.tree)
        #cur.tree.cpr<-pcrp_tree(cur.tree, 2)
        #mu_hat <- compute_mu_hat(cur.tree.cpr, 500)
        #print(mu_hat)
    }


}

apply.crp.method("cognate_ie_compatible.tree.rooted.tree")











