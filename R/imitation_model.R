#' Imitation Model For Nucleotide Sequences
#'
#' @description The function creates a list of nucleotide sequences based on the
#'     first generated sequence with specified parameters. Imitation model
#'     consists of two parts. *generate_positions()* function returns nucleotide
#'     positions that will be regenerated, while the *generate_sequences()*
#'     function performs generation of sequences. As a result, the obtained
#'     sequences will have some similarity to each other, leading them to be
#'     clustered together.
#'
#' @param number_of_sequences Integer, number of sequences in each cluster
#' @param sequences_length Integer, the length of each sequence
#' @param percent_of_new_positions Double, percentage of nucleotides that will
#'     be regenerated in each sequence
#' @param nucleotides Character, alphabet of symbols from which sequences are
#      formed
#' @param number_of_clusters Integer, number of clusters into which the data
#'     will be divided
#'
#' @return The function returns a list of nucleotide sequences.
#' @export
#'
#' @examples
#' # imitation_model(6, 15, 0.4, c("A", "T", "C", "G"), 3)
#' # imitation_model(100, 200, 0.9, c("A", "T", "C", "G"), 12)

imitation_model <- function(number_of_sequences, sequences_length, percent_of_new_positions, nucleotides, number_of_clusters){
  replicate(number_of_clusters, {
    sequences <- vector("list", number_of_sequences)
    sequence <- replicate(sequences_length, sample(nucleotides, 1, replace = TRUE))
    sequences[[1]] <- paste(sequence, collapse = "")

    positions <- generate_positions(number_of_sequences, sequences_length, percent_of_new_positions)
    sequences <- generate_sequences(sequences, sequences_length, positions, nucleotides)
  })
}

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
