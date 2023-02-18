fourier_series <- function(k, coef, x) {
  # Compute the values of the Fourier series at the given points
  y <- apply(matrix(1:k, ncol = k), 2, function(j) {
    cos_term <- cos(2 * pi * (j - 1) * x)
    sin_term <- sin(2 * pi * (j - 1) * x)
    coef[j] * cos_term + coef[j + k] * sin_term
  })

  # Return the values of the Fourier series
  return(rowSums(y))
}



# KL-represent the covariance function
KL_fourier <- funciton(t,s, k=50, nu = -2, scale = 1) {

  seq_vec <- (1:k)^(nu) * scale

  fourier_f <-

}


# Define the length of the sequence
N <- 10

# Initialize an empty vector to store the sequence

sum(seq_vec)



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

