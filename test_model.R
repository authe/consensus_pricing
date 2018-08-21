# estimate model with simulated data
source('EstimateModel2_Lagged_pubpriv.R')
source('SimulateData2_Lagged_pubpriv.R')

S <- 15
TT <- 250
n.sim <- 1
ord <- 2
paras.sim <- c(0.7, -1.4, 0.1, 0.01, 0.03)

l <- 1e-10
h <- 1

est.sim <- c()

for (seed in 1:n.sim){
  print(seed)
  yt.sim <- SimulateDataLagged(S = S, TT = TT, rho = paras.sim[1], theta.bar = paras.sim[2], sig.u = paras.sim[3], sig.v = paras.sim[4], sig.n = paras.sim[5], ord = ord, seed = seed)
  
  # initial value for parameters
  paras.0 <- c(0.5, -1, 0.1, 0.1, 0.1)
  
  # initialise prior mean and covariance of state vector
  a0 <- c(rep(0, ord+1), rep(0,S*ord), 0)
  P0 <- var(yt.sim[1,]) * diag(((1+ord) + S*ord + 1))
  
  # P0 <- diag(S) %x% sp$SIG
  # P0 <- cbind( matrix(0, ncol=(ord+1), nrow=(S*ord)), P0)
  # P0 <- rbind( cbind(sp$PP, matrix(0, nrow=(ord+1), ncol=(S*ord))), P0)
  # P0 <- cbind(P0, rep(0,(ord+1) + (S*ord)))
  # P0 <- rbind(P0, c(rep(0,(ord+1) + (S*ord)), sp$PP[2,2]))
  
  est.out <- estimate.model(paras0 = paras.0, yt = yt.sim, a0 = a0, P0 = P0, ord = 2, l = l, h = h)
  est.sim <- rbind(est.sim, c(est.out$par, est.out$LL, est.out$converge))
}

colnames(est.sim) <- c('rho', 'theta', 'sig.u', 'sig.v', 'sig.n', 'LL', 'con')