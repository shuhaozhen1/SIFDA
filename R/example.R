

covmatrix <- diag(c(1,1,1,1,1)); covmatrix[4,5] <- 0.5 ; covmatrix[5,4] <- 0.5

t_points <- seq(0,1,0.1)

coef_list <- list(
  beta1 = function(t){0.5*cos(2*pi*(t-0.5))},
  beta2 = function(t){cos(2*pi*t)},
  beta3 = function(t){sin(2*pi*t)},
  beta4 = function(t){4*t*(1-t)},
  # beta5 = function(t){0.4*(t-0.6)^2}
  beta5 = function(t){2*(t-0.6)^2}
)

cov1 <- function(x,y) {
  k <- 100

  fourier <- 0
  for (i in 1:k) {
    fourier <- fourier + i^(-2) * (sin(pi*i *x) + cos(pi*i * x)) *  (sin(pi*i *y) + cos(pi*i * y))
  }

  return(fourier)

}


cov_list <- list(
  x1 = function(t,s){2^(t/2+s/2-abs(t-s)/4)},
  x2 = function(t,s){exp(t/2+s/2-abs(t-s)/4)},
  x3 = function(t,s){exp(t/2+s/2-abs(t-s)/4)},
  x4 = function(t,s){exp(-abs(t-s)/4)},
  x5 = function(t,s){exp(-abs(t-s)/4)}
)

mean_list <- rep(list(function(x){0*x}),5)


#VCMdata <- rp_VCM_generate(50,5, coef_list, mean_list, cov_list)

#VCM_inference(VCMdata,t_points = t_points)

#
#
# est_VCM(VCMdata, t=t_points , h=0.1)-VCM_inference(VCMdata,t_points = t_points, h=0.1)$betahat
