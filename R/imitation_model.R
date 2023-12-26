generate_positions <- function(number_of_sequences, sequences_length, percent_of_new_positions){
  num_samples <- as.integer(percent_of_new_positions * sequences_length)
  positions <- replicate(num_samples, sample(1:sequences_length, 1))
  unlist(positions)
}

generate_sequences <- function(sequences, sequences_length, positions, nucleotides) {
  base_sequence <- strsplit(sequences[[1]], "")[[1]]
  modified_sequences <- lapply(sequences[-1], function(seq) {
    new_sequence <- base_sequence
    new_sequence[positions] <- sample(nucleotides, sum(length(positions)), replace = TRUE)
    paste(new_sequence, collapse = "")
  })
  c(sequences[1], modified_sequences)
}

imitation_model <- function(number_of_sequences, sequences_length, percent_of_new_positions, nucleotides, number_of_clusters){
  replicate(number_of_clusters, {
    sequences <- vector("list", number_of_sequences)
    sequence <- replicate(sequences_length, sample(nucleotides, 1, replace = TRUE))
    sequences[[1]] <- paste(sequence, collapse = "")

    positions <- generate_positions(number_of_sequences, sequences_length, percent_of_new_positions)
    sequences <- generate_sequences(sequences, sequences_length, positions, nucleotides)
  })
}
