SP.model <- function(rho, theta.bar, sig.u, sig.v, sig.n, ord, S, sig.z=0, tol=1e-15){

  #########################  CREATE STATE SPACE MODEL FOR FKF PACKAGE ########################################
  
  # (i) bring model into state space form for FKF package Kalman filter
  # s_t = F s_(t-1) + H eps_t, eps_t ~ N(0,I)
  # y_t = Z d_t + G gamma_t , gamma_t ~ N(0,I)
  # here s_t = (theta_t^(0), ..., theta_t^(ord), x_(1,t)^(1),...,x_(1,t)^(ord), ..., x_(S,t)^(ord), theta_{t-1}^(1))
  # where x_(i,t)^(l) = theta_(i,t)^(l) - theta_t^(l) and ord >= 1
  # y_t = (p_t, theta_(1,t)^(1), theta_(2,t)^(1),...,theta_(S,t)^(1))
  
  source("createStateSpaceMat2_Lagged_pubpriv.R")
  
  aux <- StateSpaceMatLag(ord = ord, rho = rho, sig.u = sig.u, sig.v = sig.v, sig.n = sig.n, tol = tol)
  M <- aux$M
  N <- aux$N
  Q <- aux$Q
  V <- aux$V
  PP <- aux$PP  # used in initial condition P0 for fkf
  SIG <- V %*% t(V) # used in initial condition P0 for fkf
  
  # transition equation: alpha_t = dt + Tt alpha_{t-1} + Ht eta_t
  # drift
  dt <- matrix(0, nrow = (ord + 1 + ((ord)*S) + 1))
  
  # transition matrix
  Tt <- diag(S) %x% Q     # %x% is Kronecker product
  Tt <- cbind(Tt, matrix(0, ncol=1, nrow=(S*ord)))
  Tt <- cbind(matrix(0, ncol=(ord+1), nrow=(S*ord)), Tt)
  Tt <- rbind(Tt, c(0,1, rep(0, (ord-1) + ord*S + 1) ) )  # unit vector to select theta_{t-1}^(1)
  Tt <- rbind( cbind(M, matrix(0, ncol = ((S*ord) + 1), nrow=(ord+1) )), Tt)
  
  # shock variance matrix HHt = E(Ht eta eta' Ht') = Ht Ht'
  Ht <- diag(S) %x% V
  Ht <- cbind(matrix(0, ncol=2, nrow=(S*ord)), Ht)
  Ht <- rbind(Ht, rep(0, 2 + S))  # row of zeros for theta_{t-1}^(1)
  Ht <- rbind( cbind(N, matrix(0, ncol=S, nrow=(ord+1))), Ht)
  
  HHt <- Ht %*% t(Ht)
  
  # observation equation: y_t = ct + Zt alpha_t + Gt epsilon_t
  # drift
  ct <- matrix(theta.bar, nrow = (S+1))
  
  # observation matrix
  Zt <- matrix(0, nrow=(1+S), ncol=(ord+1 + (S*ord) + 1))
  Zt[1,ord+1 + (ord*S) + 1] <- 1  
  Zt[2:(S+1),2] <- 1
  for (s in 1:S){
    Zt[(1+s),( ord+1 + (s-1)*(ord) + 1)] <- 1
  }
  
  # shock variance matrix GGt = E(Gt e e' Gt') = Gt Gt'
  Gt = diag(sig.z, S)
  Gt = rbind(rep(0,S), Gt)
  
  GGt <- Gt %*% t(Gt)
  #GGt <- diag(0,N+1)
  
  return(list(dt = dt, Tt = Tt, HHt = HHt, ct = ct, Zt = Zt, GGt = GGt, PP = PP, SIG = SIG))
}