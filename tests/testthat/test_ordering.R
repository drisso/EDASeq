context("Test same ordering of rows after normalization if no rownames.")

counts <- matrix(rnbinom(n=1e3, size=2,
                         mu=rnorm(n=1e3, mean=100, sd=10)),
                 nrow=100, ncol=10)
gc <- runif(n=100, min=0.3, max=0.8)

test_that("ordering is correct without rownames", {
  # normalize a matrix without rownames (and strip any name given by the function)
  normCounts_noNames <- unname(withinLaneNormalization(unname(counts), gc, which="full"))
  # normalize a matrix with rownames (and strip any name given by the function)
  rownames(counts) <- paste0("gene", 1:nrow(counts))
  normCounts_withNames <- unname(withinLaneNormalization(counts, gc, which="full"))

  expect_equal(normCounts_noNames, normCounts_withNames)
}
)
