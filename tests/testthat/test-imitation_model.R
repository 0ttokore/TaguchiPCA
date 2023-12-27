test_that("generate_positions returns correct data type", {
  number_of_sequences = 50
  sequences_length = 100
  percent_of_new_positions = 0.4
  positions <- generate_positions(number_of_sequences, sequences_length, percent_of_new_positions)
  expect_type(positions, "integer")
  expect_length(positions, sequences_length * percent_of_new_positions)
})

test_that("generate_sequences returns correct data type and length", {
  sequences <- replicate(10, paste(sample(c("A", "T", "C", "G"), 50, replace = TRUE), collapse = ""))
  positions <- sample(1:50, 20)
  modified_sequences <- generate_sequences(sequences, 50, positions, c("A", "T", "C", "G"))
  expect_type(modified_sequences, "list")
  expect_length(modified_sequences, 10)
})
