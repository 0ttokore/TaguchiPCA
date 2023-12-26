-*
TaguchiClustering <- function(sequences, number_of_clusters){
  dist_matrix <- stringdist::stringdistmatrix(as.character(sequences), as.character(sequences), method = 'lv')^2

  row_means <- rowMeans(dist_matrix)
  col_means <- colMeans(dist_matrix)
  grand_mean <- mean(dist_matrix)

  X_centered <- sweep(sweep(dist_matrix, 1, row_means, "-"), 2, col_means, "-")
  X_centered <- (X_centered + grand_mean)*-1/2

  result <- eigen(X_centered)

  vectors <- result$vectors[,1:10]
  values <- matrix(ncol = 10, nrow = 1)
  for(i in 1:10){
    values[,i] <- result$values[[i]][1]
  }

  PC <- matrix(ncol = 10, nrow = nrow(X_centered))
  for(i in 1:10){
    PC[,i] <- sqrt(values[,i]) * vectors[,i]
  }

  mDist <- dist(scale(X_centered),method = "euclidian")
  Wardlustering <- hclust(mDist,method = "ward.D2")

  plot(Wardlustering)
  rect.hclust(Wardlustering, k = number_of_clusters)

  groups <- cutree(Wardlustering, k = number_of_clusters)
  plot(PC, main = "PCoA components:", col = groups,lwd = 2, xlab = 'PC1', ylab = 'PC2')
  grid()

  pairs(PC[,1:10], col = groups, lwd = 2, main = "Clustered data in the new PCoA feature space:")
  grid()

  K <- 2
  Sb <- matrix(0, dim(X_centered)[1], dim(X_centered)[1])

  for (i in 1:K) {
    v <- vectors[,i]
    v_outer <- tcrossprod(v, v)
    SB <- values[i] * v_outer
    Sb <- Sb + SB
  }

  Znew <- X_centered %*% Sb

  pairs(Znew[,1:10],col = groups, main = "Clustered data in extended feature space:")
  plot(Znew[,1:2], main = "Clustered data in extended feature space", col = groups,lwd = 1, xlab = 'PC1\'', ylab = 'PC2\'')
  grid()
}
