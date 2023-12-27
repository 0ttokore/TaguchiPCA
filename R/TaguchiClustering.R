#' Taguchi Clustering Method
#'
#' @description The function takes a list of sequences and a distance metric,
#'     based on which the given sequences will be vectorized. Then the function
#'     computes and returns the *PC_num* number of principal coordinates and also
#'     performs a transition to a new feature space.
#'
#' @param sequences List, non-vectorized nucleotide sequences
#' @param number_of_clusters Integer, number of clusters which you want to divide
#'     your sequences. The default is 4
#' @param method Method for distance calculation. The default is "lv", but you can
#'     also use another method: *'osa'* – optimal alignment, *'dl'* –
#'     Damerau-Levenshtein, *'qgram'* – q gram, *'cosine'* – cosine, *'jw'* – Jaro,
#'     *'lcs'* – greatest common substrings
#' @param PC_num Integer, number of principal coordinates u want to take from the
#'     function
#' @param clusters_plot Bool, displays a plot on the first two principal coordinates
#' @param new_clusters_plot Bool, displays a plot on the first two principal
#'     coordinates in new feature space
#'
#' @return List of *P_num* principal coordinates, *P_num* coordinates in the new
#'     feature space
#' @export
#' @importFrom graphics grid
#' @importFrom graphics pairs
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats rect.hclust
#' @importFrom stringdist stringdistmatrix
#'
#' @examples
#' TaguchiClustering(imitation_model(), number_of_clusters = 8)
#' TaguchiClustering(imitation_model(), method = 'qgram', new_clusters_plot = TRUE)

TaguchiClustering <- function(sequences, method = 'lv', number_of_clusters = 4, PC_num = 10, clusters_plot = FALSE, new_clusters_plot = FALSE){
  dist_matrix <- stringdist::stringdistmatrix(as.character(sequences), as.character(sequences), method = method)^2

  row_means <- rowMeans(dist_matrix)
  col_means <- colMeans(dist_matrix)
  grand_mean <- mean(dist_matrix)

  X_centered <- sweep(sweep(dist_matrix, 1, row_means, "-"), 2, col_means, "-")
  X_centered <- (X_centered + grand_mean)*-1/2

  result <- eigen(X_centered)

  vectors <- result$vectors[,1:PC_num]
  values <- matrix(ncol = PC_num, nrow = 1)
  for(i in 1:PC_num){
    values[,i] <- result$values[[i]][1]
  }

  PC <- matrix(ncol = PC_num, nrow = nrow(X_centered))
  for(i in 1:PC_num){
    PC[,i] <- sqrt(values[,i]) * vectors[,i]
  }

  if(clusters_plot == TRUE || new_clusters_plot == TRUE){
    mDist <- dist(scale(X_centered),method = "euclidian")
    Wardlustering <- hclust(mDist,method = "ward.D2")
    groups <- cutree(Wardlustering, k = number_of_clusters)
  }

  if(clusters_plot == TRUE){
    plot(Wardlustering)
    rect.hclust(Wardlustering, k = number_of_clusters)

    plot(PC, main = "PCoA components:", col = groups,lwd = 2, xlab = 'PC1', ylab = 'PC2')
    grid()
  }

  K <- 2
  Sb <- matrix(0, dim(X_centered)[1], dim(X_centered)[1])

  for (i in 1:K) {
    v <- vectors[,i]
    v_outer <- tcrossprod(v, v)
    SB <- values[i] * v_outer
    Sb <- Sb + SB
  }

  Znew <- X_centered %*% Sb

  if(new_clusters_plot == TRUE){
    pairs(Znew[,1:10],col = groups, main = "Clustered data in extended feature space:")

    plot(Znew[,1:2], main = "Clustered data in extended feature space", col = groups, lwd = 1, xlab = 'PC1', ylab = 'PC2')
    grid()
  }

  return(list(PC = PC, Znew = Znew))
}
