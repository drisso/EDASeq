context("Test plot functions.")

set.seed(42)

mat <- matrix(data=rpois(1000, lambda=10), ncol=10)
rownames(mat) <- paste("gene", 1:nrow(mat), sep="")
colnames(mat) <- paste("sample", 1:ncol(mat), sep="")

es <- newSeqExpressionSet(mat)

test_that("plotPCA works", {
  ks <- 2:5

  ## matrix
  expect_silent(lapply(ks, function(k) plotPCA(mat, k=k)))
  expect_silent(lapply(ks, function(k) plotPCA(mat, k=k, labels=FALSE)))
  expect_silent(lapply(ks, function(k) plotPCA(mat, k=k)))
  expect_silent(lapply(ks, function(k) plotPCA(mat, k=k, labels=FALSE, pch=20, col=1:2)))

  ## expressionset
  expect_silent(lapply(ks, function(k) plotPCA(es, k=k)))
  expect_silent(lapply(ks, function(k) plotPCA(es, k=k, labels=FALSE)))
  expect_silent(lapply(ks, function(k) plotPCA(es, k=k, labels=FALSE, pch=20, col=1:2)))

})

test_that("plotRLE works", {

  ## matrix
  rle <- plotRLE(mat)
  expect_equal(dim(mat), dim(rle))

  ## expressionset
  rle <- plotRLE(es)
  expect_equal(dim(mat), dim(rle))

})

