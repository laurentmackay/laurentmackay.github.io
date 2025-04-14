compute_AUCs<-function(df, C_thresh, dose=1){
  df$t_thresh <- sapply(with(df, log(dose/(C_thresh*V))/k), function(x) max(x,0))
  df$auc_infinity <- with(df, dose/(V*k))
  df$auc_thresh <- with(df, auc_infinity*(1.0-exp(-k*t_thresh)) )
  df
}


# log_concentration_sample<-function(k,V){
#   n<-10
#   t_sample<-0.5+cumsum(runif(n, 0.02,0.03))
#   data.frame(
#     time=t_sample,
#     log_concentration=-k*t_sample-log(V)+rnorm(n,0,0.2),
#   )
# }

concentration_sample<-function(k,V, t_sample){
  n<-length(t_sample)
  data.frame(
    time=t_sample,
    concentration=exp(-k*t_sample+rnorm(n,0,0.2))/V
  )
}