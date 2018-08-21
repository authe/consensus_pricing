estimate.model <- function(paras0, yt, a0, P0, ord, l, h, opt.algo="Nelder-Mead"){

  # ML estimation of model given data and initial conditions
  # parameters to estimate: rho (persistence fundamental), sig.u (fundamental shock), sig.e (public noise), sig.n (private noise), sig.z (measurement error)
  # paras0: vector of initial conditions c(rho,theta.bar,sig.u,sig.e,sig.n,sig.z [optional]) to be used in optimization
  # a0, P0: prior mean and covariance matrix for state vector
  # additional parameters:
  #   l, h : upper and lower bounds for sig.u, sig.e, sig.n, sig.z
  #   ord: order of model
  #   opt.algo : optimisation algorithm to be used by optimr ("Nelder_Mead" by default)
  #
  # returns: parameter estimates, log-likelihood value, and convergence indicator of optimr
  
  # load libraries and auxiliary functions
  library(tictoc)
  library(optimr)
  source('LogLike2_Lagged_pubpriv.R')

  ####################################### objective function passed to 'optim'  ###############################

  objective <- function(paras, yt, ord, l, h, tol=1e-15) {
    # ord: beliefs order used
    # l: lower bound for standard deviation
    # h : upper bound for standard deviations
    
    rho <- 1/(1+exp(-paras[1]))
    theta.bar <- paras[2]
    sig.u <- (h + l*exp(-paras[3]))/(1 + exp(-paras[3]))
    sig.v <- (h + l*exp(-paras[4]))/(1 + exp(-paras[4]))
    sig.n <- (h + l*exp(-paras[5]))/(1 + exp(-paras[5]))
    
    par <- c(rho, theta.bar, sig.u, sig.v, sig.n)
    
    LL.val <- LL.model(paras = par, yt = yt, a0 = a0, P0 = P0, ord = ord, tol = tol)
    
    return(LL.val)
  }
  
  # objective function evaluated at data yt 
  f <- function(x){
    val <- objective(paras=x, yt=yt, ord = ord, l = l, h = h, tol = 1e-15)
    return(val)
  }

  ###################################  estimation  ##########################################################
  
  # transform initial conditions for optimisation
  
  trans.rho.0 <- log(paras0[1]/(1-paras0[1]))
  theta.bar.0 <- paras0[2]
  log.sig.u.0 <- log((paras0[3]-l)/(h-paras0[3]))
  log.sig.v.0 <- log((paras0[4]-l)/(h-paras0[4]))
  log.sig.n.0 <- log((paras0[5]-l)/(h-paras0[5]))
  paras0.trans <- c(trans.rho.0, theta.bar.0, log.sig.u.0, log.sig.v.0, log.sig.n.0)

  # optimisation using opt.algo
  #fit <- optimr(par = paras0.trans, fn = f, method=opt.algo, hessian = TRUE, control = list(maxit = 10000))
  fit <- optim(par = paras0.trans, fn = f, method=opt.algo, hessian = TRUE, control = list(maxit = 10000))
  
  est <- fit$par
  estimates <- c(rho=1/(1+exp(-est[1])), theta.bar=est[2], sig.u=(h + l*exp(-est[3]))/(1 + exp(-est[3])), sig.v=(h + l*exp(-est[4]))/(1 + exp(-est[4])), sig.n=(h + l*exp(-est[5]))/(1 + exp(-est[5])))
  
  out.data <- list(par=estimates, LL=fit$value, converge=fit$convergence)
  
  return(out.data)
  
}