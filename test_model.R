# estimate model with simulated data
source('EstimateModel2_Lagged_pubpriv.R')
source('SimulateData2_Lagged_pubpriv.R')

set.seed(1)

n.test <- 2
n.sim <- 10

S <- 5
TT <- 500
ord <- 2
l <- 1e-10
h <- 1

#paras.sim <- c(0.95, -1.4, 0.1, 0.15, 0.1)
paras.sim <- matrix(rnorm(n.test*5), ncol=5, nrow=n.test)
paras.sim[,1] <- 1/(1+exp(-paras.sim[,1]))
paras.sim[,3] <- (h + l*exp(-paras.sim[,3]))/(1 + exp(-paras.sim[,3]))/10
paras.sim[,4] <- (h + l*exp(-paras.sim[,4]))/(1 + exp(-paras.sim[,4]))/10
paras.sim[,5] <- (h + l*exp(-paras.sim[,5]))/(1 + exp(-paras.sim[,5]))/10

est.list <- list()
length(est.list) <- n.test

for (n in 1:n.test){
  est.sim <- c()
  
  for (seed in 1:n.sim){
    print(seed)
    yt.sim <- SimulateDataLagged(S = S, TT = TT, rho = paras.sim[n,1], theta.bar = paras.sim[n,2], sig.u = paras.sim[n,3], sig.v = paras.sim[n,4], sig.n = paras.sim[n,5], ord = ord, seed = seed)
    
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
    aux <- c(rep(rho.0 * sig2.theta.0, (1+ord)), rep(0,ord*S))
    P0 <- cbind(P0, aux)
    P0 <- rbind(P0, c(aux,sig2.theta.0))
    
    est.out <- estimate.model(paras0 = paras.0, yt = yt.sim, a0 = a0, P0 = P0, ord = ord, l = l, h = h)
    est.sim <- rbind(est.sim, c(est.out$par, est.out$LL, est.out$converge))
    colnames(est.sim) <- c('rho', 'theta', 'sig.u', 'sig.v', 'sig.n', 'LL', 'con')
    print(est.sim)
  }
  est.list[[n]] <- est.sim 
}