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
KL_fourier_cov <- function(t,s) {
  k<-100
  seq_vec <- (1:k)^(-2)

  fourier_f <- sapply(1:k, function(i,x,y,c){
    c[i] * (sin(2*pi*i *x) + cos(2*pi*i * x)) *  (sin(2*pi*i *y) + cos(2*pi*i * y))
  },x=t,y=s,c=seq_vec)

  return(sum(fourier_f))

}

KL_fourier_cov(0.4,0.61)


outer(seq(0,1,0.1), KL_fourier_cov)

outer(x,x, function(x,y){KL_fourier_cov(x,y,k=50, nu = -2, scale = 1)})


# Define the covariance function with fixed arguments
cov_func <- function(x, y, a, b) {
  exp(-a * (x - y)^2) + b
}

# Generate a vector of points
x <- seq(0, 1, length.out = 5)

# Compute the covariance matrix using outer
outer(x, x, function(x, y) cov_func(x, y, a = 1, b = 0.5))




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

