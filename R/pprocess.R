#' Function to generate values of a random process at specified time points
#'
#' @param time_points: a vector of time points at which to evaluate the random process
#' @param p: number of columns in the output, representing different realizations of the random process,
#'  to be interpreted as dimensions. Default is 5.
#' @param meanf: a scalar value representing the mean of the random process
#' @param covariancef: a covariance function that maps two time points to their covariance
#' @param num_basis: the number of basis functions to use for approximation. Default is 1000.
#' @param distribution: the distribution to use for generating weights for the basis functions. Default is "normal".
#'   Accepted values are "normal", "uniform", and "exponential".
#
#' @return A matrix with `length(time_points)` rows and `p+1` columns, the first column records the time points and
#' each other column represents a realization of the random process at the specified time points.
#' @export
generate_random_process_values <- function(time_points, meanf=function(x){0},
                                           covariancef, num_basis=1000, distribution = "normal") {

  # Calculate the covariance matrix between each pair of time points
  covariance_matrix <- outer(time_points, time_points, covariancef)

  # Decompose the covariance matrix into eigenvectors and eigenvalues
  eigen <- eigen(covariance_matrix)

  # Determine the number of basis functions to use for approximation
  # If the requested number of basis functions is greater than the number of eigenvectors, use all eigenvectors instead
  num_basis <- min(num_basis, ncol(eigen$vectors))

  # Select the first `num_basis` eigenvectors as the basis functions
  basis_functions <- eigen$vectors[, 1:num_basis]
  eigen_values <- eigen$values[1:num_basis]

  # Generate weights for the basis functions
  # The distribution for generating the weights can be specified using the `distribution` argument
  if (distribution == "normal") {
    weights <- rnorm(num_basis)
  } else if (distribution == "uniform") {
    weights <- runif(num_basis)
  } else if (distribution == "exponential") {
    weights <- rexp(num_basis)-1
  } else {
    stop("Invalid distribution argument. Choose 'normal', 'uniform', or 'exponential'.")
  }

  # Scale the weights by the square root of the eigenvalues
  weights <-  t(sqrt(diag(eigen_values)) %*% weights)

  # Calculate the values of the random process at each time point
  # by taking the inner product between the weights and the basis functions
  process_values_matrix <- weights %*% t(basis_functions)
  process_values_matrix <- t(process_values_matrix) + meanf(time_points)

  # Return the matrix of random process values
  return(cbind(time_points, process_values_matrix))
}

pprocess_tp <- function(time_points, mean_list, cov_list, n_b= 1000, dis_vector= NULL, covmatrix=NULL){

  p <- length(mean_list)

  if(is.null(dis_vector)){
    dis_vector <- rep('normal',p)
  }

  if(is.null(covmatrix)){
    covmatrix <- diag(rep(1),p)
  }

  x <- matrix(0, nrow = length(time_points), ncol = p)
  for( i in 1:p) {
    x[,i] <- generate_random_process_values(time_points=time_points, meanf=mean_list[[i]],
                                   covariancef=cov_list[[i]],
                                   num_basis=n_b, distribution = dis_vector[i])[,2]
  }

  trans_x <- t(covmatrix %*% t(x))

  return(cbind(t=time_points, xp = trans_x))
}



### FDA
rp_generate <- function(n,m,mean_list,cov_list,dis_vector=NULL,covmatrix=NULL,sig=0.1,depend=F, domain=c(0,1)){
  generate_time_points <- function(n, m, domain = c(0,1)) {
    # Generate a list to store the time points for each element
    time_points_list <- list()

    for (i in 1:n) {
      # Sample the number of time points for the i-th element from (m-1, m, m+1)
      m_i <- sample(c(m-1, m, m+1), 1)

      # Generate m_i time points uniformly distributed from 0 to 1
      time_points_i <- sort(runif(m_i, min = domain[1], max = domain[2]))

      # Add the time points for the i-th element to the list
      time_points_list[[i]] <- time_points_i
    }

    # Return the list of time points
    return(time_points_list)
  }

  t_list <- generate_time_points(n=n,m=m, domain = domain)

  p_list <- lapply(t_list, pprocess_tp, mean_list=mean_list, cov_list=cov_list,dis_vector=dis_vector,
                         covmatrix=covmatrix)


  if(depend == F){
    p_error <- lapply(p_list, function(x){
      x[,-1] <- x[,-1] + rnorm(n=nrow(x), sd= sig)
      return(x)
      })
  } else {
    p_error <- lapply(p_list, function(x){
      x[,-1] <- x[,-1] + rnorm(1, sd= sig) * x[,2]/length(mean_list)
      return(x)
    })
  }

  return(p_error)
}

### VCM generation
rp_VCM_generate <- function(n,m,coef_list, mean_list,cov_list,dis_vector=NULL,
                            covmatrix=NULL,sig=0.1,depend=F, domain=c(0,1)){
  generate_time_points <- function(n, m, domain = c(0,1)) {
    # Generate a list to store the time points for each element
    time_points_list <- list()

    for (i in 1:n) {
      # Sample the number of time points for the i-th element from (m-1, m, m+1)
      m_i <- sample(c(m-1, m, m+1), 1)

      # Generate m_i time points uniformly distributed from 0 to 1
      time_points_i <- sort(runif(m_i, min = domain[1], max = domain[2]))

      # Add the time points for the i-th element to the list
      time_points_list[[i]] <- time_points_i
    }

    # Return the list of time points
    return(time_points_list)
  }

  t_list <- generate_time_points(n=n,m=m, domain = domain)

  p_list <- lapply(t_list, pprocess_tp, mean_list=mean_list, cov_list=cov_list,dis_vector=dis_vector,
                   covmatrix=covmatrix)

  for(i in 1:n) {
    beta <- sapply(coef_list, function(x){x(p_list[[i]][,1])})
    yi <-  rowSums( p_list[[i]][,-1] * beta )
    p_list[[i]] <- cbind(p_list[[i]],yi)
  }

  if(depend == F){
    p_error <- lapply(p_list, function(x){
      x[,ncol(x)] <- x[,ncol(x)] + rnorm(n=nrow(x), sd= sig)
      return(x)
    })
  } else {
    p_error <- lapply(p_list, function(x){
      x[,ncol(x)] <- x[,ncol(x)] + rnorm(1, sd= sig) * x[,ncol(x)]/length(mean_list)
      return(x)
    })
  }

  return(p_error)
}


