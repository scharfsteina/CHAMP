find_min <- function(gamma_matrix) {
  n <- dim(gamma_matrix)[1]
  mins <- c()
  partitions <- c()
  i = 1 # location in matrix
  while (i < n) {
    curr_gamma <- min(gamma_matrix[i,(i+1):n], na.rm=T)
    mins  <- append(mins, curr_gamma)
    loc <- which(gamma_matrix==curr_gamma, arr.ind=T)
    partitions <- append(partitions, i)
    i = loc[1]
  }
  return(list("mins" = mins, "partitions" = partitions))
}