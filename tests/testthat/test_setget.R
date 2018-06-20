context("Test plot functions.")

set.seed(42)

mat <- matrix(data=rpois(1000, lambda=10), ncol=10)
rownames(mat) <- paste("gene", 1:nrow(mat), sep="")
colnames(mat) <- paste("sample", 1:ncol(mat), sep="")

es <- newSeqExpressionSet(mat)

test_that("setters and getters works", {

  counts(es) <- mat
  expect_equal(counts(es), mat)

  normCounts(es) <- mat
  expect_equal(normCounts(es), mat)
})

test_that("between-lane normalization works", {

  es <- newSeqExpressionSet(mat)
  es2 <- newSeqExpressionSet(mat)

  # full
  norm_mat <- betweenLaneNormalization(mat, which="full")
  normCounts(es2) <- norm_mat
  es3 <- betweenLaneNormalization(es, which="full")

  expect_equal(es2, es3)

  # median
  norm_mat <- betweenLaneNormalization(mat, which="median")
  normCounts(es2) <- norm_mat
  es3 <- betweenLaneNormalization(es, which="median")

  expect_equivalent(es2, es3)

  # upper
  norm_mat <- betweenLaneNormalization(mat, which="upper")
  normCounts(es2) <- norm_mat
  es3 <- betweenLaneNormalization(es, which="upper")

  expect_equivalent(es2, es3)

})

test_that("within-lane normalization works", {

  gc <- rnorm(100)
  es <- newSeqExpressionSet(mat, featureData = data.frame(gc=gc, row.names=rownames(mat)))
  es2 <- es

  # full
  norm_mat <- withinLaneNormalization(mat, gc, which="full")
  normCounts(es2) <- norm_mat
  es3 <- withinLaneNormalization(es, "gc", which="full")

  expect_equal(es2, es3)

  # median
  norm_mat <- withinLaneNormalization(mat, gc, which="median")
  normCounts(es2) <- norm_mat
  es3 <- withinLaneNormalization(es, "gc", which="median")

  expect_equivalent(es2, es3)

  # upper
  norm_mat <- withinLaneNormalization(mat, gc, which="upper")
  normCounts(es2) <- norm_mat
  es3 <- withinLaneNormalization(es, "gc", which="upper")

  expect_equivalent(es2, es3)

  # loess
  norm_mat <- withinLaneNormalization(mat, gc, which="loess")
  normCounts(es2) <- norm_mat
  es3 <- withinLaneNormalization(es, "gc", which="loess")

  expect_equivalent(es2, es3)
})
