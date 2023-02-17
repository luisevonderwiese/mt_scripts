library(phangorn)

cog.tree=read.tree("IE2011_Cognates_rel_ANNOT.nwk")
tr = read.tree("../output/project_pavlos/pars_results/pars_current_opt_tree_large_step10000000.nwk")
d=treedist(tr, cog.tree, check.labels=TRUE)
print(paste("***** d: ", d[1]))
