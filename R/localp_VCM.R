est_t_VCM <- function(data,h=0.1,t,d=1){
  n<- length(data)
  p <- ncol(data[[1]])-2
  mis <- sapply(data, function(x){nrow(x)})

  totaldata  <- Reduce(rbind, data)

  totalt <- totaldata[,1]

  totalxp <- totaldata[,2:(p+1)]

  h_design <- diag(rep(1,p)) %x% diag(h^(0:d))

  totalxp_design <- matrix(0,nrow = nrow(totaldata), ncol= p*(d+1))
  for(k in seq(1,by=d+1, length.out=p)){
    for(i in 0:d)
      totalxp_design[,k+i] <- totalxp[,(k-1)/(d+1)+1] * (totalt-t)^i
  }

  K_design <- diag(Epa_K((totalt-t)/h)/h/(rep(mis,mis)))

  y_design <- totaldata[,p+2]

  est <- (solve(t(totalxp_design) %*% K_design %*% totalxp_design) %*%
    t(totalxp_design) %*% K_design %*% y_design)[seq(1,by=d+1, length.out=p),]

  return(est)
}


#' @title Estimate the coefficients of a VCM model
#' @description Given a list of data, estimate the coefficients of a VCM model at given time points
#'
#' @param data A list of data frames, each of which represents the observed data for one element.
#' @param t_points A vector of time points at which to estimate the VCM coefficients.
#' @param h The bandwidth for the Epanechnikov kernel. Default is 0.1.
#' @param d The degree of the polynomial for the local linear regression. Default is 1.
#'
#' @return A matrix of estimated coefficients, with each row representing the estimated coefficients for one element.
#'
#' @examples
#' data_list <- rp_VCM_generate(n=3,m=10,coef_list=list(f1,f2),mean_list=list(mean1,mean2),
#'                     cov_list=list(cov1,cov2), sig=0.1)
#' t_points <- seq(0, 1, by = 0.1)
#' est <- est_VCM(data_list, t_points)
#' est
#'
#' @export
est_VCM <- function(data,t_points, h=0.1,d=1,h1=NULL) {
  est <- sapply(t_points, est_t_VCM,h=h,d=d, data=data)
  return(est)
}


### center data

center_VCM <- function(data,h=0.1,d=1){
  totaldata <- Reduce(rbind, data)
  totalxp <- totaldata[,2:(ncol(totaldata)-1)]

  mis <- sapply(data, function(x){nrow(x)})

  estbeta <- t(est_VCM(data=data, t_points = totaldata[,1]))

  epsilon <- totaldata[,ncol(totaldata)] - rowSums(estbeta * totalxp)

  totaldata[,ncol(totaldata)] <- epsilon

  centerdata <- split_matrix_by_rows(totaldata, mis)

  return(centerdata)
}



### individual for VCM
localp_t_VCM <- function(data_i,data,h=0.1,t,d=1){
  n<- length(data)
  p <- ncol(data[[1]])-2
  mis <- sapply(data, function(x){nrow(x)})

  totaldata  <- Reduce(rbind, data)

  totalt <- totaldata[,1]

  totalxp <- totaldata[,2:(p+1)]

  totaly <- totaldata[,p+2]

  totalxp_design <- matrix(0,nrow = nrow(totaldata), ncol= p*(d+1))
  for(k in seq(1,by=d+1, length.out=p)){
    for(i in 0:d)
    totalxp_design[,k+i] <- totalxp[,(k-1)/(d+1)+1] * (totalt-t)^i
  }

  K_design <- diag(Epa_K((totalt-t)/h)/h/(rep(mis,mis)))

  h_design <- diag(rep(1,p)) %x% diag(h^(0:d))

  denominator <- solve( (1/n) * solve(h_design) %*%
    t(totalxp_design) %*% K_design %*% totalxp_design %*% solve(h_design)
    )

  xt <- data_i[,1]
  xp <- data_i[,2:(p+1)]
  y <- data_i[, p+2]

  xp_design <- matrix(0,nrow = length(y), ncol= p*(d+1))
  for(k in seq(1,by=d+1, length.out=p)){
      for(i in 0:d)
        xp_design[,k+i] <- xp[,(k-1)/(d+1)+1] * (xt-t)^i
  }


  if(is.null(h1)){
    h1 <- 0.6*h
  }
  h1_design <- diag(rep(1,p)) %x% diag((h1)^(0:d))

  k_d <- diag(Epa_K((xt-t)/h1)/h1/length(y))


  xi <- (solve(h1_design)%*% denominator %*% solve(h1_design)%*%
             t(xp_design) %*% k_d %*% y)[seq(1,by=d+1, length.out=p),]

  return(xi)
}

localp_ts_VCM <- function(data_i,data,h=0.1,h1=NULL,t_points,d=1){
  xi_t <- sapply(t_points, localp_t_VCM, data_i=data_i, data= data, h=h, h1=h1, d=d)
  return(xi_t)
}

localp_VCM_i <- function(data, h=0.1, h1=NULL, t_points, d=1) {
  xi <- lapply(X=data, localp_ts_VCM, data= data, h=h, h1=h1, t_points=t_points, d=d)
  return(xi)
}


#' @title Estimate and infer variance component model
#'
#' @description This function estimates the variance component model using local polynomial regression and performs inference for the coefficients.
#'
#' @param data A list of data frames, each containing data for one cluster.
#' @param bstime Number of bootstrap samples for confidence interval estimation.
#' @param h The bandwidth for the local polynomial regression. If NULL, the bandwidth will be chosen using cross-validation.
#' @param t_points A vector of time points to evaluate the estimate.
#' @param d The degree of the local polynomial regression.
#' @param alpha The significance level for the confidence interval.
#'
#' @return A list containing the estimate of the coefficients (betahat), the lower and upper bounds of the confidence interval (low and up), the average bandwidth (band), and the estimated standard deviation (sd).
#' \item{betahat}{An estimate of the coefficient functions at the specified time points.}
#' \item{low}{The lower bound of the confidence band.}
#' \item{up}{The upper bound of the confidence band.}
#' \item{band}{The average width of the confidence band.}
#' \item{sd}{The standard deviation of the estimated function.}
#' @examples
#' data(VCMdata)
#' VCM_inference(data = VCMdata, bstime = 100, h = NULL, t_points = seq(0, 10, by = 0.1), d = 1, alpha = 0.05)
#'
#' @importFrom MASS mvrnorm
#' @export

VCM_inference <- function(data, bstime=3000, h=NULL, h1=NULL, t_points, d=1, alpha = 0.05){

  n <- length(data)
  p <- ncol(data[[1]])-2

  if(is.null(h)){
    totalt <- Reduce(rbind,data)[,1]
    h <- bw.ucv(totalt)
  }

  betahat <- est_VCM(data = data, t_points = t_points, h=h, d=d)

  centerdata <- center_VCM(data = data, h=h,d=d)

  xis <- localp_VCM_i(data=centerdata,h=h, t_points = t_points, d=d)

  xis_v <- lapply(xis, function(x){c(t(x))})

  xis2_ts <- lapply(xis_v, function(x){x %*% t(x)})

  est_cov <- Reduce('+',xis2_ts)/n

  est_sd <- sqrt(diag(est_cov))

  standard <- est_cov / (est_sd %*% t(est_sd))


  bs <- abs(MASS::mvrnorm(n=bstime, mu=rep(0,length(t_points)*p), Sigma = standard))

  q_hat <- quantile(apply(bs, 1, max), 1-alpha)

  est_sd <- t(matrix(est_sd, ncol=p))

  scb_low <- betahat - q_hat * est_sd / sqrt(n)
  scb_up <- betahat + q_hat * est_sd / sqrt(n)

  return(list(betahat= betahat, low=scb_low, up=scb_up,
              band= mean(scb_up-scb_low), sd=est_sd))
}



