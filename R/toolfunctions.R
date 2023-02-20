

# Define the Epanechnikov kernel function
Epa_K <- function(x) {
  # Use sapply to apply the function to each element of x
  y <- sapply(x, function(xi) {
    # Check if xi is within [-c, c]
    if (abs(xi) <= 1) {
      # If xi is within [-c, c], return the Epanechnikov kernel value
      return(0.75 * (1 - xi^2))
    } else {
      # If xi is outside [-c, c], return 0
      return(0)
    }
  })
  return(y)
}




split_matrix_by_rows <- function(mat, lengths) {
  cumulative_lengths <- c(0, cumsum(lengths))
  submatrices <- lapply(1:(length(cumulative_lengths) - 1), function(i) mat[(cumulative_lengths[i] + 1):cumulative_lengths[i + 1], ])
  submatrices
}

