context("Test GC-content and length.")

test_that("getGeneLengthAndGCContent works", {
  gcc <- getGeneLengthAndGCContent(id=c("ENSG00000012048", "ENSG00000139618"), org="hsa", mode = "biomart")

  expect_equal(nrow(gcc), 2)
}
)
