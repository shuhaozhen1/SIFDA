

cv_h_VCM <- function(data, h, k=5, d=1){
  n <- length(data)
  leaveout <- seq(1,n,by=k)

  p <- ncol(data[[1]])-2

  mse <- 0

  for (i in leaveout){

    new_data_list <- data[-c(i:(i+k-1))]

    leavedata <- Reduce(rbind,data[c(i:(i+k-1))])

    leavet <- leavedata[,1]

    leavey <- leavedata[,p+2]

    new_est <- est_VCM(data = new_data_list, t_points = leavet,h=h, d=d)

    new_y <- rowSums(leavedata[,2:(p+1)] * t(new_est))

    mse <- mse + sum((leavey - new_y)^2)
  }

  return(mse)
}

cv_VCM <- function(data,h_seq,k=5,d=1){
  mses <- sapply(h_seq,cv_h_VCM,data=data,k=k,d=d)
  h_cv <- h_seq[which.min(mses)]
  return(h_cv)
}



true_h_VCM <- function(data, h, d=1,coef_list) {
  t <- seq(0,1,0.01)
  est <- est_VCM(data = data, t_points = t,h=h, d=d)
  truebeta <- sapply(coef_list, function(f,x){f(x)},x=t)

  mse <- sum((est - t(truebeta))^2)/length(t)

  return(mse)
}






