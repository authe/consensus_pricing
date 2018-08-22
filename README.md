# consensus_pricing

## List of file with description

* **createStateSpaceMat2_Lagged_pubpriv.R**: Defines the function _StateSpaceMatLag(ord, rho, sig.u, sig.v, sig.n)_ that creates matrices for the transition and observation equation of the submitter's state space system.
* **FKFModel2_Lagged_pubpriv.R** : Using *createStateSpaceMat2_Lagged_pubpriv.R* as an input, defines function _SP.model(rho, theta.bar, sig.u, sig.v, sig.n, ord, S, sig.z=0)_ that creates matrices for econometrician's state space system. Outputs all necessary matrices to compute the loglikehood with the Kalman Filter implementation of the FKF package.
* **LogLike2_Lagged_pubpriv.R** : Defines function _LL.model(paras, yt, a0, P0, ord, tol=1e-15)_ that takes parameter vector (rho,theta,sig.u,sig.v,sig.n), data yt, and prior mean a0 and covariance P0 of state vector as input and returns loglikelihood of data and posterior covariance matrix of state vector.
* **SimulateData2_Lagged_pubpriv** : Simulates consensus price and individual submissions given parameters (rho,theta,sig.u,sig.v,sig.n).
* **EstimateModel2_Lagged_pubpriv.R** : Estimates parameters of model via maximum likelihood using _optim_ given data yt and prior mean a0 and covariance P0.
* **test_model.R** : Estimates parameter of models using simulated data.