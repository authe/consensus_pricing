StateSpaceMatLag <- function(ord, rho, sig.u, sig.v, sig.n, tol=1e-15){

  # function to create matrices for state-state state space system with beliefs of order ord
  # output: M_k, N_k, PP_k, KK_k 
  # transition equation: theta_t = M_k theta_{t-1} + N_k w_t with w_t ~ N(0, I_2)
  # observation equation: y_t = D_{k,1} theta_t + D_{k,2} theta_{t-1} + R_w w_t + R_n eta_{j,t} with eta_{j,t} = N(0,1)
  # where theta_t = (theta_t^(0), theta_t^(1), ..., theta_t^(ord)) and ord >= 1
  # y_(1,t) = theta_t + sig.n eta_{j,t} + sig.v v_t and y_(2,t) = theta_{t-1}^(1)
  # P_k = lim_t -> infty Var(theta_t) : stationary state covariance matrix
  # KK_k : stationary Kalman gain
  
    
  ##########################  SET INITIAL CONDITIONS   #########################################
  
  # initialise matrices for iteration of state space system
  # theta_t = A theta_(t-1) + C (w_t eta_{j,t})'
  # y_t = D1 theta_t + D2 theta_{t-1} + R (w_t eta_{j,t})'
  # (formulas used for Kalman filter derive from Nimark (2015) "A Low Dimensional Kalman Filter", Economic Letters)
  
  # initialise consensus price signal at p_t = theta_{t-1}^(0) + sig.v w_{2,t}
  
  # initial conditions for state equation: theta_t^(0) = rho theta_{t-1}^(0) + sig.u w_{1,t} 
  A <- matrix(rho, nrow=1)
  C <- matrix( c(sig.u, 0, 0), nrow=1)
  
  # initial conditions for observation equation
  D1 <- matrix(c(1,0), nrow=2)   
  D2 <- matrix(c(0,1), nrow=2)
  R.w <- matrix( c(0,0,sig.v,sig.v), nrow=2)
  R.n <- matrix( c(sig.n, 0), nrow=2 )
  R <- cbind(R.w, R.n)
  
  # starting value covariance matrix
  P <- matrix( sig.u^2/(1- rho^2), nrow=1)
  
  ########################## ITERATE to get transition matrices with next higher order of beliefs  ######################
  
  for (k in 1:(ord)){
    
    # update Kalman gain KK and stationary state covariance matrix P
    
    #tol <- 10^(-15)
    count.max <- 100
    
    test.convergence <- TRUE
    count <- 0
    while (test.convergence & count <= count.max ){
      KK <- ( A %*% P %*% t(D1 %*% A + D2) + C %*% t(C) %*% t(D1) + C %*% t(R)) %*% solve( (D1 %*% A + D2) %*% P %*% t(D1 %*% A + D2) + (D1 %*% C + R) %*% t(D1 %*% C + R))
      PP <- A %*% P %*% t(A) + C %*% t(C) - KK %*% ( (D1 %*% A + D2) %*% P %*% t(D1 %*% A + D2) + (D1 %*% C + R) %*% t(D1 %*% C + R)) %*% t(KK)
      
      # matrix distance between PP and P
      dist <- max(abs(PP - P))
      test.convergence <- (dist > tol)
      
      P <- PP
      count <- count + 1
      
    }
    
    # update state system
    
    # record old values (output for individual transition matrices)
    A.old <- A
    D1.old <- D1
    D2.old <- D2
    KK.old <- KK
    
    # update A, C, and D
    
    R.w <- matrix( c(0,0,sig.v,0), nrow=2)
    
    MA <- rbind( c(rho, rep(0,k)), matrix(0, nrow=k, ncol=k+1) )
    MB <- rbind( rep(0,k), KK %*% (D1 %*% A + D2) )
    MB <- cbind( MB, rep(0,k+1))
    MC <- rbind( rep(0, k), A - KK %*% (D1 %*% A + D2) )
    MC <- cbind( rep(0, k+1), MC)
    A <- MA + MB + MC
    
    C <- KK %*% (D1 %*% matrix(C[,1:2], ncol=2) +  R.w)
    C <- rbind( sig.u * c(1,0), C)
    C <- cbind(C, rep(0,k+1))
    
    D1 <- cbind(D1, c(0,0))
    D2 <- rbind(rep(0, k+1), c(0, 1, rep(0,k-1)))
    
    P <- cbind(P, rep(0,k))
    P <- rbind(P, rep(0,k+1))
    P[k+1,k+1] = P[k,k]
    
  }
  
  # update Kalman gain KK and stationary state covariance matrix P
  #tol <- 10^(-15)
  count.max <- 100
  
  test.convergence <- TRUE
  count <- 0
  while (test.convergence & count <= count.max ){
    KK <- ( A %*% P %*% t(D1 %*% A + D2) + C %*% t(C) %*% t(D1) + C %*% t(R)) %*% solve( (D1 %*% A + D2) %*% P %*% t(D1 %*% A + D2) + (D1 %*% C + R) %*% t(D1 %*% C + R))
    PP <- A %*% P %*% t(A) + C %*% t(C) - KK %*% ( (D1 %*% A + D2) %*% P %*% t(D1 %*% A + D2) + (D1 %*% C + R) %*% t(D1 %*% C + R)) %*% t(KK)
    
    # matrix distance between PP and P
    dist <- max(abs(PP - P))
    test.convergence <- (dist > tol)
    
    P <- PP
    count <- count + 1
    
  }
  
  # calculate matrices for individual deviations from theta_t
  # x_(j,t) = Q_k x_(j,t-1) + V_k eta_(j,t)
  
  Q <- A.old - KK.old %*% (D1.old %*% A.old + D2.old)
  V <- matrix(sig.n * KK.old[,1], ncol=1)
  
  
  ###############################   return results    ####################################################
  
  KalMat <- list()
  KalMat$"M" <- A
  KalMat$"N" <- matrix(C[,1:2], ncol=2)
  KalMat$"PP" <- PP
  KalMat$"KK" <- KK.old
  KalMat$"M.ind" <- A.old
  KalMat$"D1" <- D1.old
  KalMat$"D2" <- D2.old
  KalMat$"Q" <- Q
  KalMat$"V" <- V
  
  return(KalMat)
  
}