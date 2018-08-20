# estimate model with simulated data

source('EstimateModel2_Lagged_pubpriv.R')
source('SimulateData2_Lagged_pubpriv.R')

S <- 15
TT <- 250
ord <- 2
paras.sim <- c(0.95, -1.3, 0.2, 0.01, 0.04)

yt.sim <- SimulateDataLagged(S = S, TT = TT, rho = paras.sim[1], theta.bar = paras.sim[2], sig.u = paras.sim[3], sig.v = paras.sim[4], sig.n = paras.sim[5], ord = ord)

l <- 1e-10
h <- 1

paras.0 <- c(0.54, -1.45, 0.124, 0.108, 0.104)

est.out <- estimate.model(paras0 = paras.0, yt = yt.sim, ord = 2, l = l, h = h)