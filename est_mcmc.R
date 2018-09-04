library('mcmc')
library('tictoc')

source('LogLike2_demeaned_Lagged_pubpriv.R')
source('SimulateData2_Lagged_pubpriv.R')

ll.mcmc <- function(par, yt, a0, P0, ord, l, h){
  # transform parameters
  paras <- rep(0,4)
  paras[1] <- ifelse(par[1] > 0, 1/(1+exp(-par[1])), exp(par[1])/(1+exp(par[1])))
  paras[2] <- ifelse(par[2] > 0, (h + l*exp(-par[2]))/(1 + exp(-par[2])), (l + h*exp(par[2]))/(1+exp(par[2])))
  paras[3] <- ifelse(par[3] > 0, (h + l*exp(-par[3]))/(1 + exp(-par[3])), (l + h*exp(par[3]))/(1+exp(par[3])))
  paras[4] <- ifelse(par[4] > 0, (h + l*exp(-par[4]))/(1 + exp(-par[4])), (l + h*exp(par[4]))/(1+exp(par[4])))
  
  ll.out <- LL.model.demeaned(paras = paras, yt = yt, a0 = a0, P0 = P0, ord = ord)
  
  # calculate Jacobian for transformed parameters
  # (see https://stats.stackexchange.com/questions/330402/what-is-an-example-of-a-transformation-on-a-posterior-distribution-such-that-the )
  log.jac <- rep(0,4)
  log.jac[1] <- ifelse(par[1] > 0, -log1p(exp(-par[1]))-log1p(exp(par[1])), par[1] - 2*log1p(exp(par[1])))
  log.jac[2] <- ifelse(par[2] > 0, log(h-l) -log1p(exp(-par[2]))-log1p(exp(par[2])), log(h-l) + par[2] - 2*log1p(exp(par[2])))
  log.jac[3] <- ifelse(par[3] > 0, log(h-l) -log1p(exp(-par[3]))-log1p(exp(par[3])), log(h-l) + par[3] - 2*log1p(exp(par[3])))
  log.jac[4] <- ifelse(par[4] > 0, log(h-l) -log1p(exp(-par[4]))-log1p(exp(par[4])), log(h-l) + par[4] - 2*log1p(exp(par[4])))
  
  ll <- ll.out$LL + sum(log.jac) - 3*log(h-l)
      
  return(ll)
}


# simulated data

S <- 20
TT <- 500
ord <- 2
seed <- 1

l <- 0
h <- 5

paras.sim <- c(0.9, -1.4, 0.1, 0.02, 0.1)

yt.sim <- SimulateDataLagged(S = S, TT = TT, rho = paras.sim[1], theta.bar = paras.sim[2], sig.u = paras.sim[3], sig.v = paras.sim[4], sig.n = paras.sim[5], ord = ord, seed = seed)

# initial value for parameters
consIV <- yt.sim[1,]
est.AR1 <- arima(consIV, c(1,0,0))
rho.0 <- min(abs(est.AR1$coef[1]),0.99)
sig.u.0 <- sqrt(est.AR1$sigma2)
aux <- sweep(yt.sim[2:(S+1),1:(TT-1)], 2, yt.sim[1,2:TT])
sig.n.0 <- mean(apply(aux, 2, sd, na.rm=TRUE))
sig.v.0 <- sig.n.0

# initialise prior mean and covariance of state vector
a0 <- c(rep(0, ord+1), rep(0,S*ord), 0)
P0 <- diag(S) %x% matrix(sig.n.0^2, nrow=ord, ncol=ord)
P0 <- rbind(matrix(0, nrow=(1+ord), ncol=(S*ord)), P0)
sig2.theta.0 <- sig.u.0^2 /(1 - rho.0^2)
aux <- matrix(sig2.theta.0, nrow=(1+ord), ncol=(1+ord))
aux <- rbind(aux, matrix(0, nrow=(S*ord), ncol=(1+ord)))
P0 <- cbind(aux, P0)
aux <- c(rep(rho.0 * sig2.theta.0, (1+ord)), rep(0,ord*S))
P0 <- cbind(P0, aux)
P0 <- rbind(P0, c(aux,sig2.theta.0))

yt.sim.demeaned <- yt.sim - mean(yt.sim)

#par.0 <- c(0.7, 0.1, 0.2, 0.1)
par.0 <- c(rho.0, sig.u.0, sig.v.0, sig.n.0)
par.0 <- log(par.0/(1- par.0))

# Rprof('prof.txt')
# tic()
# ll.mcmc(par = rnorm(4), yt = yt.sim, a0 = a0, P0 = P0, ord = ord, l=0, h=1)
# toc()
# Rprof(NULL)
# summaryRprof('prof.txt')

tic()
out.mcmc <- metrop(ll.mcmc, initial = par.0, nbatch = 3e3, scale = 0.03, yt = yt.sim.demeaned, a0 = a0, P0 = P0, ord = ord, l=l, h=h)
toc()