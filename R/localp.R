library(MASS)

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

  denominator <- solve(t(totalxp_design) %*% K_design %*% totalxp_design)

    xt <- data_i[,1]
    xp <- data_i[,2:(p+1)]
    y <- data_i[, p+2]

    xp_design <- matrix(0,nrow = length(y), ncol= p*(d+1))
    for(k in seq(1,by=d+1, length.out=p)){
      for(i in 0:d)
        xp_design[,k+i] <- xp[,(k-1)/(d+1)+1] * (xt-t)^i
    }

    k_d <- diag(Epa_K((xt-t)/h)/h/length(y))

    xi <- n * (denominator %*% t(xp_design) %*% k_d %*% y)[seq(1,by=d+1, length.out=p),]

  return(xi)
}

localp_ts_VCM <- function(data_i,data,h=0.1,t_points,d=1){
  xi_t <- sapply(t_points, localp_t_VCM, data_i=data_i, data= data, h=h, d=d)
  return(xi_t)
}

localp_VCM_i <- function(data, h=0.1, t_points, d=1) {
  xi <- lapply(data, localp_ts_VCM, data= data, h=h, t_points=t_points, d=d)
  return(xi)
}

VCM_inference <- function(data, bstime=1000, h=NULL, t_points, d=1, alpha = 0.05){

  n <- length(data)
  p <- ncol(data[[1]])-2

  if(is.null(h)){
    totalt <- Reduce(rbind,data)[,1]
    h <- 2* bw.bcv(totalt)
  }

  xis <- localp_VCM_i(data,h= h, t_points = t_points, d=d)

  betahat <- Reduce('+',xis)/n

  xis_v <- lapply(xis, function(x){c(t(x))})

  xis2_ts <- lapply(xis_v, function(x){x %*% t(x)})
  est_cov <- Reduce('+',xis2_ts)/n - c(t(betahat)) %*% t(c(t(betahat)))

  est_sd <- sqrt(diag(est_cov))

  standard <- est_cov / (est_sd %*% t(est_sd))

  bs <- abs(mvrnorm(n=bstime, mu=rep(0,length(t_points)*p), Sigma = standard))

  q_hat <- quantile(apply(bs, 1, max), 1-alpha)

  est_sd <- t(matrix(est_sd, ncol=p))

  scb_low <- betahat - q_hat * est_sd / sqrt(n)
  scb_up <- betahat + q_hat * est_sd / sqrt(n)

  return(list(betahat= betahat, scb_low=scb_low, scb_up=scb_up, band= mean(scb_up-scb_low)/2))
}


