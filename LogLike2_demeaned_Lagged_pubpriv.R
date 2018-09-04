# returns log-likelihood of model computed by FKF package
LL.model.demeaned <- function(paras, yt, a0, P0, ord, tol=1e-15){
  
  # a0 and P0 are prior mean and covariance matrix for state vector
 
  library(FKF)
  source("FKFModel2_demeaned_Lagged_pubpriv.R")
  
  S <- dim(yt)[1] - 1
  TT <- dim(yt)[2]
  
  rho <- paras[1]
  sig.u <- paras[2]
  sig.v <- paras[3]
  sig.n <- paras[4]
  
  sp <- SP.model.demeaned(rho = rho, sig.u = sig.u, sig.v = sig.v, sig.n = sig.n, ord = ord, tol = tol , S = S)
  
  ans <- fkf(a0 = a0, P0 = P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
  
  out.data <- list(LL=ans$logLik, Ptt=ans$Ptt[,,TT])
  
  return(out.data)
   
}