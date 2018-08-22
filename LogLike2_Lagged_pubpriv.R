# returns log-likelihood of model computed by FKF package
LL.model <- function(paras, yt, a0, P0, ord, tol=1e-15){
  
  # a0 and P0 are prior mean and covariance matrix for state vector
 
  library(FKF)
  source("FKFModel2_Lagged_pubpriv.R")
  
  S <- dim(yt)[1] - 1
  TT <- dim(yt)[2]
  
  rho <- paras[1]
  theta.bar <- paras[2]
  sig.u <- paras[3]
  sig.v <- paras[4]
  sig.n <- paras[5]
  
  sp <- SP.model(rho = rho, theta.bar = theta.bar, sig.u = sig.u, sig.v = sig.v, sig.n = sig.n, ord = ord, tol = tol , S = S)
  
  ans <- fkf(a0 = a0, P0 = P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
  
  out.data <- list(LL=ans$logLik, Ptt=ans$Ptt[,,TT])
  
  return(out.data)
   
}