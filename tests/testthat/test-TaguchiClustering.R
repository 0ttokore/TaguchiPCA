test_that("TaguchiClustering returns correct data type, objects and data dimensions", {
  fake_sequences <- replicate(100, paste(sample(c("A", "T", "C", "G"), 50, replace = TRUE), collapse = ""))
  result <- TaguchiClustering(fake_sequences)

  expect_type(result, "list")

  expect_true("PC" %in% names(result))
  expect_true("Znew" %in% names(result))

  expect_equal(nrow(result$PC), 100)
  expect_equal(ncol(result$PC), 10)
  expect_equal(nrow(result$Znew), 100)
  expect_equal(ncol(result$Znew), 100)
})
