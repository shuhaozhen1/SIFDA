est_t_Mean <- function(data,h=0.1,t,d=1){
  n<- length(data)
  p <- ncol(data[[1]])-1
  mis <- sapply(data, function(x){nrow(x)})

  totaldata  <- Reduce(rbind, data)

  totalt <- totaldata[,1]

  totalyp <- totaldata[,2:(p+1)]

  totalt_design <- matrix(0,nrow = nrow(totaldata), ncol=(d+1))

  for(k in 1:(d+1)){
    totalt_design[,k] <- (totalt-t)^(k-1)
  }

  K_design <- diag(Epa_K((totalt-t)/h)/h/(rep(mis,mis)))

  est <- (solve(t(totalt_design) %*% K_design %*% totalt_design) %*%
            t(totalt_design) %*% K_design %*% totalyp)[1,]

  return(est)
}

#' @title Estimate the mean function at give time points
#'
#' @description This function estimates the mean function at given time poitns using local polynomial regression estimator with an Epanechnikov kernel.
#'
#' @param data A list of matrices, where each matrix contains the data for a single subject.
#' @param h The bandwidth parameter for the kernel function.
#' @param t_points The time points at which the mean function are estimated.
#' @param d The degree of the polynomial used in the local polynomial regression.
#'
#' @return A matrix representing the estimated mean function at the given time points.
#'
#' @examples
#' data <- rp_generate(n=10, m=20, mean_list=list(function(x) (x+1)^2, function(x) (x-1)^2),
#'             cov_list=list(function(x,y) exp(-abs(x-y)), function(x,y) exp(-2*abs(x-y))),
#'             sig=0.1)
#' est_Mean(data, t_points = seq(0,1,0.1))
#'
#' @export
est_Mean <- function(data,t_points, h=0.1,d=1) {
  est <- sapply(t_points, est_t_Mean,h=h,d=d, data=data)
  return(est)
}

center_Mean <- function(data,h=0.1,d=1){
  totaldata <- Reduce(rbind, data)

  totalt <- totaldata[,1]

  mis <- sapply(data, function(x){nrow(x)})

  esty <- t(est_Mean(data=data, t_points = totaldata[,1]))

  epsilon <- totaldata[,-1] - esty

  totaldata[,-1] <- epsilon

  centerdata <- split_matrix_by_rows(totaldata, mis)

  return(centerdata)
}


### individual for Mean
localp_t_Mean <- function(data_i,data,h=0.1,h1=NULL,t,d=1){
  n<- length(data)
  p <- ncol(data[[1]])-1
  mis <- sapply(data, function(x){nrow(x)})

  totaldata  <- Reduce(rbind, data)

  totalt <- totaldata[,1]

  totalyp <- totaldata[,2:(p+1)]

  h_design <- diag(h^(0:d))

  totalt_design <- matrix(0,nrow = nrow(totaldata), ncol=(d+1))

  for(k in 1:(d+1)){
    totalt_design[,k] <- (totalt-t)^(k-1)
  }

  K_design <- diag(Epa_K((totalt-t)/h)/h/(rep(mis,mis)))

  denominator <- solve( (1/n) * solve(h_design) %*%
                          t(totalt_design) %*% K_design %*% totalt_design %*% solve(h_design) )

  ti <- data_i[,1]
  ypi <- data_i[,2:(p+1)]

  ti_design <- matrix(0,nrow = nrow(data_i), ncol=(d+1))

  for(k in 1:(d+1)){
    ti_design[,k] <- (ti-t)^(k-1)
  }

  if(is.null(h1)){
    h1 <- 0.5*h
  }

  ki <- diag(Epa_K((ti-t)/h1)/h1/length(ti))

  h1_design <- diag(h1^(0:d))

  ihat <- (solve(h1_design) %*% denominator %*% solve(h1_design)%*%
           t(ti_design) %*% ki %*% ypi)[1,]

  return(ihat)
}


localp_ts_Mean <- function(data_i,data,h=0.1,h1=NULL,t_points,d=1){
  xi_t <- sapply(t_points, localp_t_Mean, data_i=data_i, data= data, h=h, h1=h1, d=d)
  return(xi_t)
}


localp_Mean_i <- function(data, h=0.1,h1=NULL, t_points, d=1) {
  xi <- lapply(X=data, localp_ts_Mean, data= data, h=h, h1=h1, t_points=t_points, d=d)
  return(xi)
}


#' @title Mean Inference for Local Polynomial Regression
#'
#' @description This function performs mean inference for local polynomial regression.
#'
#' @param data A list of data frames, where each data frame contains the time series data for one subject. The first column of each data frame should contain the time points and the remaining columns should contain the observations at each time point.
#' @param bstime The number of bootstrap replicates.
#' @param h The bandwidth parameter for the local polynomial regression. If not provided, the bandwidth will be selected using the biweight cross-validation method.
#' @param t_points The time points at which to estimate the mean function.
#' @param d The degree of the local polynomial regression.
#' @param alpha The level of the confidence interval.
#'
#' @return A list containing the following components:
#' \item{muhat}{An estimate of the mean function at the specified time points.}
#' \item{low}{The lower bound of the confidence interval.}
#' \item{up}{The upper bound of the confidence interval.}
#' \item{band}{The average width of the confidence interval.}
#' \item{sd}{The standard deviation of the estimated mean function.}
#'
#' @examples
#' data <- rp_generate(n=10, m=20, mean_list=list(function(x) (x+1)^2, function(x) (x-1)^2),
#'             cov_list=list(function(x,y) exp(-abs(x-y)), function(x,y) exp(-2*abs(x-y))),
#'             sig=0.1)
#'
#' Mean_inference(data = data, t_points = seq(0,1,0.1))
#'
#' @importFrom MASS mvrnorm
#'
#' @export
Mean_inference <- function(data, bstime=1000, h=NULL, h1=NULL, t_points, d=1, alpha = 0.05){

  n <- length(data)
  p <- ncol(data[[1]])-1

  if(is.null(h)){
    totalt <- Reduce(rbind,data)[,1]
    h <- 1.2 * bw.bcv(totalt)
  }

  muhat <- est_Mean(data = data, t_points = t_points, h=h, d=d)

  centerdata <- center_Mean(data = data, h=h,d=d)

  if(is.null(h1)){
    h1 <- 0.5 * h
  }

  xis <- localp_Mean_i(data=centerdata,h= h1, t_points = t_points, d=d)

  xis_v <- lapply(xis, function(x){c(t(x))})

  xis2_ts <- lapply(xis_v, function(x){x %*% t(x)})

  est_cov <- Reduce('+',xis2_ts)/n

  est_sd <- sqrt(diag(est_cov))

  standard <- est_cov / (est_sd %*% t(est_sd))

  bs <- abs(MASS::mvrnorm(n=bstime, mu=rep(0,length(t_points)*p), Sigma = standard))

  q_hat <- quantile(apply(bs, 1, max), 1-alpha)

  est_sd <- t(matrix(est_sd, ncol=p))

  scb_low <- muhat - q_hat * est_sd / sqrt(n)
  scb_up <- muhat + q_hat * est_sd / sqrt(n)

  return(list(muhat= muhat, low=scb_low, up=scb_up,
              band= mean(scb_up-scb_low), sd=est_sd))
}



