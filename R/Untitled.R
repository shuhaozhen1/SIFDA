# # Load the ggplot2 library
# library(ggplot2)
#
# # Create some example data
# VCMdata <- rp_VCM_generate(200,5, coef_list, mean_list, cov_list)
#
# df <- VCM_inference(VCMdata,t_points = t_points,h=0.35)
#
#
# df <- data.frame(t= t_points, hat = df$betahat[1,], low = df$scb_low[1,],
#                  up= df$scb_up[1,], true = coef_list[[1]](t_points))
#
#
# ggplot(df, aes(x = t)) +
#   geom_line(aes(y=true),color ='blue' )+
#   geom_line(aes(y=hat),color ='red' )+
#   geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.2)
#
#
# VCMdata <- rp_VCM_generate(50,5, coef_list, mean_list, cov_list)
#
# hatn <- function(n,m,h){
#   VCMdata <- rp_VCM_generate(n,m, coef_list, mean_list, cov_list)
#   hat <- est_VCM(VCMdata,t_points = t_points, h=h)
#   return(hat)
# }
#
# a <- replicate(100, hatn(15,50,0.15),simplify = F)
#
# sqrt(20 * (Reduce('+',lapply(a, function(x){x^2}))/100 - (Reduce('+',a)/100 )^2))
#
# xis <- localp_VCM_i(VCMdata, t_points = t_points, d=1)
#
# apply(Reduce(rbind,VCMdata),2,max)
#
# xis2 <- lapply(xis,function(x){x^2})
#
# Reduce('+',lapply(xis, function(x){x^2}))/50 - (Reduce('+',xis)/50 )^2
#
# xis <- localp_VCM_i(VCMdata, h=0.15, t_points = t_points)
# xis2 <- lapply(xis,function(x){x^2})
# Reduce('+',lapply(xis, function(x){x^2}))/50 - (Reduce('+',xis)/50 )^2
# hat <- Reduce('+',localp_VCM_i(VCMdata, h=0.1, t_points = t_points))/50
#
# wildbsonetime <- function(data,t_points,h,d=1){
#   n<- length(data)
#   resample <- sample(1:n,n,replace = T)
#   newdata <- data[resample]
#   newest <- est_VCM(data=newdata,t_points = t_points,h=h,d=d)
#   return(newest)
# }
#
# a <- replicate(100,wildbsonetime(VCMdata,t_points,h=0.15,d=1),simplify = F)
# sqrt(50 * (Reduce('+',lapply(a, function(x){x^2}))/100 - (Reduce('+',a)/100 )^2))
