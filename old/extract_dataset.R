a <- read.table("Datasets without geo.xlsx - Full Dataset.csv", h=TRUE, sep=",", quote="\"")
colnames(a)[colnames(a) == "Ptg.E"] <- "PtgE"
glos <- read.table('glosses.txt')[,1]

## check if this holds
glos[!(glos %in% colnames(a))]

data <- a[,glos]
data <- t(data)
colnames(data) <- a[glos %in% colnames(a),3]

write.table(data, file="46glossesfull.csv", sep="\t", row.names=T, col.names=T, quote=F)
