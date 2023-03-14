library(CRPTree)


apply.crp.method <- function(tree.name) {
    data = read.table("46glossesfull.csv", header=TRUE, row.names=1, colClasses='character')
    taxa_names = row.names(data)
    datamat = as.matrix(data)
    tab <- apply(datamat, 2, function(x){table(factor(x, levels=c(0,1)))})
    indsToUse <- apply(tab, 2, function(x){x[1] > 1 && x[2] > 1})
    datamat <- datamat[,indsToUse]
    tree<-read.tree(file=tree.name)
    tree$tip.label[tree$tip.label == "Ptg-E"] = "PtgE"
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
        class(cur.tree) <- 'phylo'
        cur.tree<-process_tree(cur.tree)
        cur.tree.cpr<-pcrp_tree(cur.tree, 2)
        mu_hat <- compute_mu_hat(cur.tree.cpr, 500)
        print(mu_hat)
    }


}

#apply.crp.method("IE2011_Cognates_rel_ANNOT.nwk")
apply.crp.method("geo_duration.tree")




