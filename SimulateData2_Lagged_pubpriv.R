SimulateDataLagged <- function(S, TT, rho, theta.bar, sig.u, sig.v, sig.n, ord, seed=1, sig.z=0, tol=1e-15){
  
  # code to simulate submissions for S submitters and belief order ord
  
  # import function to create state space system
  library(MASS)
  source("createStateSpaceMat2_Lagged_pubpriv.R")
  
  # simulation parameters
  set.seed(seed)     # fix seed
  # S <- 20         # number of submitters
  # TT <- 100        # submission dates
  # rho <- 0.8      # persistence of fundamental theta_t^(0) = rho theta_(t-1)^(0) + u_t  
  # sig.u <- 0.003  # sd of fundamental shock u_t
  # sig.e <- 0.01   # sd of noise shock e_t to public signal p_t = theta_t^(1) + e_t
  # sig.n <- 0.01   # sd of noise shock n_(j,t) to private signal of submitter j: s_(j,t) = theta_t^(0) + n_(j,t)
  # sig.z <- 0      # sd of measurement error for individual submissions
  # ord <- 10       # order of beliefs hierachy theta_t = (theta_t^(0),theta_t^(1),...,theta_t^(ord))
  
  # get matrices for state space system
  SSMat <- StateSpaceMatLag(ord=ord, rho=rho, sig.u=sig.u, sig.v=sig.v, sig.n=sig.n, tol=tol)
  
  #############################  SIMULATION  ##########################################################
  
  ######## simulate time series for mean beliefs
  
  # draw time series of aggregate shocks w
  w <- mvrnorm(TT+1, rep(0,2), diag(2))
  # draw initial condition using prior N(0,PP)
  theta0 <- mvrnorm(1, rep(0,ord+1), SSMat$PP)
  
  aux <- theta0
  theta <- theta0
  
  for (t in 2:(TT+1)){
    aux <- SSMat$M %*% aux + SSMat$N %*% w[t,]
    theta <- rbind(theta, t(aux))
  }
  
  
  ###### simulate time series for individual submissions
  
  # generate consensus price p_t = theta^(1)_(t-1) + eta_t
  price = c(0,theta[(1:TT),2])
  
  # generate beliefs for S submitters
  
  beliefs <- list()
 
  y0 <- mvrnorm(1, rep(0,ord+1), SSMat$PP)
  y0 <- y0[2:(ord+1)]  # initial individual beliefs
  
  for (s in 1:S){
    
    aux <- y0
    y <- aux
    
    for (t in 2:(TT+1)){
      
      priv.signal <- theta[t,1] + sig.v * w[t,2] + sig.n * rnorm(1)
      signals = matrix( c(priv.signal, price[t]), nrow=2)
      
      aux <- SSMat$M.ind %*% aux + SSMat$KK %*% ( signals - SSMat$D1 %*% SSMat$M.ind %*% aux - SSMat$D2 %*% aux) 
      y <- rbind(y, t(aux))
    }
    
    beliefs[[s]] <- y
  }
  
  # each submitter j submits his best estimate theta_(j,t)^(1)
  submissions <- c()
  for (s in 1:S){
    submissions <- cbind(submissions, beliefs[[s]][,1] + rep(theta.bar,TT+1) + sig.z * rnorm(TT+1))
  }
  
  sim.data <- t(cbind(price + theta.bar, submissions))
  
  return(sim.data)
  
}