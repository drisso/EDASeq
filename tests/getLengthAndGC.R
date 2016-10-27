library(EDASeq)
library(yeastRNASeq)

getGeneLengthAndGCContent(id=c("ENSG00000012048", "ENSG00000139618"), org="hsa", mode = "biomart")

data(geneLevelData)
getGeneLengthAndGCContent(rownames(geneLevelData)[1:3], org="sacCer3", mode="org.db")
