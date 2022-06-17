find_min1 <- function(gamma_matrix) {
  n <- dim(gamma_matrix)[1]
  mins <- c()
  partitions <- c()
  i = 1 # location in matrix
  while (i < n) {
    curr_gamma <- min(gamma_matrix[i,(i+1):n], na.rm=T)
    mins <- append(mins, curr_gamma)
    loc <- which(gamma_matrix==curr_gamma, arr.ind=T)
    loc <- loc[order(loc[,1])]
    loc <- loc[which(loc>i)]
    partitions <- append(partitions, i)
    i = loc[1]
    print(i)
  }
  return(list("mins" = mins, "partitions" = partitions))
}

find_min <- function(gamma_matrix) {
  n <- dim(gamma_matrix)[1]
  mins <- c()
  partitions <- c()
  i = 1 # location in matrix
  iprev = 1
  gamma_matrix[lower.tri(gamma_matrix)] <- NaN
  curr_gamma <- 0
  while (i <= n) {
    partitions <- append(partitions, i)
    gamma_matrix[i,(gamma_matrix[i,]<curr_gamma)] <- NaN # issue here
    if (all(is.na(gamma_matrix[i,]))) { break }
    curr_gamma <- min(gamma_matrix[i,], na.rm = T)
    mins <- append(mins, curr_gamma) # issue here
    loc <- which(gamma_matrix==curr_gamma, arr.ind=T)
    gamma_matrix[i:loc[1,2]-1,] <- NaN
    i = loc[1,2]
  }
  return(list("mins" = mins, "partitions" = partitions))
}
