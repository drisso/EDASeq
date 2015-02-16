library(EDASeq)

mat <- matrix(data=rpois(1000, lambda=10), ncol=10)
rownames(mat) <- paste("gene", 1:nrow(mat), sep="")
colnames(mat) <- paste("sample", 1:ncol(mat), sep="")

es <- newSeqExpressionSet(mat)

ks <- 2:5

## matrix
lapply(ks, function(k) plotPCA(mat, k=k))
lapply(ks, function(k) plotPCA(mat, k=k, labels=FALSE))
lapply(ks, function(k) plotPCA(mat, k=k))
lapply(ks, function(k) plotPCA(mat, k=k, labels=FALSE, pch=20, col=1:2))

## expressionset
lapply(ks, function(k) plotPCA(es, k=k))
lapply(ks, function(k) plotPCA(es, k=k, labels=FALSE))
lapply(ks, function(k) plotPCA(es, k=k, labels=FALSE, pch=20, col=1:2))
