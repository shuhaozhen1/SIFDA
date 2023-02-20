
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

#' @title Cross-validation to select bandwidth for VCM
#'
#' @description This function uses cross-validation to select the bandwidth parameter \code{h}
#' for the VCM estimation method.
#'
#' @param data a list containing matrices, each of which represents a sample of the
#' data in the VCM model. Each matrix should have a time column as the first column,
#' followed by p columns of covariates, and a final column of response values.
#' @param h_seq a sequence of candidate bandwidths to be evaluated in the
#' cross-validation.
#' @param k the number of folds in the cross-validation.
#' @param d the degree of the polynomial trend to be fitted in the VCM estimation.
#'
#' @return The value of \code{h} that minimizes the cross-validation mean squared error.
#'
#' @examples
#' h_seq <- seq(0.1,1,length=10)
#' cv_VCM(VCMdata,h_seq,k=5,d=1)
#'
#' @export
cv_VCM <- function(data,h_seq,k=5,d=1){
  mses <- sapply(h_seq,cv_h_VCM,data=data,k=k,d=d)
  h_cv <- h_seq[which.min(mses)]
  return(h_cv)
}


cv_h_Mean <- function(data, h, k=5, d=1){
  n <- length(data)
  leaveout <- seq(1,n,by=k)

  p <- ncol(data[[1]])-1

  mse <- 0

  for (i in leaveout){

    new_data_list <- data[-c(i:(i+k-1))]

    leavedata <- Reduce(rbind,data[c(i:(i+k-1))])

    leavet <- leavedata[,1]

    leaveyp <- leavedata[,2:(p+1)]

    new_est <- est_Mean(data = new_data_list, t_points = leavet,h=h, d=d)

    mse <- mse + sum((leaveyp - t(new_est))^2)
  }

  return(mse)
}

#' @title Cross-validation to select bandwidth for Mean of FDA
#'
#' @description This function uses cross-validation to select the bandwidth parameter \code{h}
#' for the mean estimation method.
#'
#' @param data a list containing matrices, each of which represents a sample of the
#' data in the FDA. Each matrix should have a time column as the first column,
#' followed by p columns response values.
#' @param h_seq a sequence of candidate bandwidths to be evaluated in the
#' cross-validation.
#' @param k the number of folds in the cross-validation.
#' @param d the degree of the polynomial trend to be fitted in the mean estimation.
#'
#' @return The value of \code{h} that minimizes the cross-validation mean squared error.
#'
#' @examples
#' h_seq <- seq(0.1,1,length=10)
#' cv_Mean(FDAdata,h_seq,k=5,d=1)
#'
#' @export
cv_Mean <- function(data,h_seq,k=5,d=1){
  mses <- sapply(h_seq,cv_h_Mean,data=data,k=k,d=d)
  h_cv <- h_seq[which.min(mses)]
  return(h_cv)
}



