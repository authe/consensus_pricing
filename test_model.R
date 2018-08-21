# estimate model with simulated data
source('EstimateModel2_Lagged_pubpriv.R')
source('SimulateData2_Lagged_pubpriv.R')

n.sim <- 10

S <- 40
TT <- 500
ord <- 2
paras.sim <- c(0.7, -1.4, 0.1, 0.01, 0.03)

l <- 1e-10
h <- 1

est.sim <- c()

for (seed in 1:n.sim){
  print(seed)
  yt.sim <- SimulateDataLagged(S = S, TT = TT, rho = paras.sim[1], theta.bar = paras.sim[2], sig.u = paras.sim[3], sig.v = paras.sim[4], sig.n = paras.sim[5], ord = ord, seed = seed)
  
  # initial value for parameters
  consIV <- yt.sim[1,]
  est.AR1 <- ar(consIV, order.max=1)
  rho.0 <- est.AR1$ar
  theta.0 <- mean(consIV)
  sig.u.0 <- sqrt(est.AR1$var.pred)
  sig.n.0 <- mean(apply(yt.sim[2:(S+1),], 2, sd, na.rm=TRUE))
  sig.v.0 <- sig.n.0
  
  paras.0 <- c(rho.0, theta.0, sig.u.0, sig.v.0, sig.n.0)
  
  # initialise prior mean and covariance of state vector
  a0 <- c(rep(0, ord+1), rep(0,S*ord), 0)
  
  #P0 <- var(yt.sim[1,]) * diag(((1+ord) + S*ord + 1))
  P0 <- diag(S) %x% matrix(sig.n.0^2, nrow=ord, ncol=ord)
  P0 <- rbind(matrix(0, nrow=(1+ord), ncol=(S*ord)), P0)
  sig2.theta.0 <- sig.u.0^2 /(1 - rho.0^2)
  aux <- matrix(sig2.theta.0, nrow=(1+ord), ncol=(1+ord))
  aux <- rbind(aux, matrix(0, nrow=(S*ord), ncol=(1+ord)))
  P0 <- cbind(aux, P0)
  aux <- c(rep(rho.0 * sig2.theta.0, 3), rep(0,ord*S))
  P0 <- cbind(P0, aux)
  P0 <- rbind(P0, c(aux,sig2.theta.0))
  
  est.out <- estimate.model(paras0 = paras.0, yt = yt.sim, a0 = a0, P0 = P0, ord = 2, l = l, h = h)
  est.sim <- rbind(est.sim, c(est.out$par, est.out$LL, est.out$converge))
}

colnames(est.sim) <- c('rho', 'theta', 'sig.u', 'sig.v', 'sig.n', 'LL', 'con')