library( doSNOW )
library( itertools )
library( doParallel )
library( foreach )
library( DFPmcmc )

dir.create("DFPSimOutput")

no_cores = 20
print(no_cores)

nrep = 50
myseed = rep(c(1:nrep), 7)
meta = c(rep(-3,nrep), rep(-2,nrep), rep(-1,nrep), rep(0,nrep), rep(1,nrep), rep(2,nrep), rep(3,nrep))
setting = data.frame(myseed,meta)

cl = makeCluster(no_cores,outfile="SimLog.txt")
registerDoSNOW(cl)
simulation = no_cores
pb <- txtProgressBar(max = simulation, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

ans <- foreach( chain_i = 1:nrow(setting), .packages = c("hDFPmcmc", "DFPmcmc", "Rcpp", "RcppArmadillo"), .options.snow = opts) %dopar%  {

  print(chain_i)
  
  N = 50 # number of subjects
  J = 5 # number of periods
  dT = 5 # col dimension of spline basis functions
  dZ = 3 # col dimension of global effect matrices
  dX = 3 # dimension of reallocation coefficients
  
  data = gendata2(mEta = setting$meta[chain_i], myseed = setting$myseed[chain_i])
  
  C = data$C
  C_ = matrix(0, nrow = N, ncol = J) # cluster indicator starts at 0
  Beta_ = array(0, dim = c(N,dT,J))
  Lambda2_ = array(1, dim = c(N,dT,J))
  Tau2_ = matrix(1, N, J)
  Nulam_ = array(1, dim = c(N,dT,J))
  Nutau_ = matrix(1, N, J)
  Gamma_ = matrix(0, N, J)
  Gamma = data$Gamma
  Eta_ = matrix(0, dX, J)
  Theta_ = matrix(0, dZ, J)
  
  D = data$D
  D_ = D*0
  
  muEta = rep(0, dX)
  muEta_ = muEta*0
  sigEta1 = 0.01*diag(dX)
  sigEta2 = 5*diag(dX)
  sigEta3 = 10*diag(dX)
  sigTheta = diag(dZ)
  
  # experiments using non-hierarchical models (with temporal dependencies)
  
  tick = Sys.time()
  set.seed(250)
  niter = 1
  DFPoutput = DFPmcmc(
    iterations = niter,
    thin = 25,
    data$YList,
    data$TList,
    data$ZList,
    data$X,
    Beta_,
    Lambda2_,
    Nulam_,
    Gamma_,
    C_,
    Tau2_,
    Nutau_,
    Eta_,
    Theta_,
    data$idstart,
    data$idend,
    SigEta = sigEta2,
    SigTheta = sigTheta,
    MuEta = muEta_,
    a = 0.1
  )
  time4 = Sys.time()-tick
  
  # experiments using non-hierarchical models (without temporal dependencies)
  
  tick = Sys.time()
  set.seed(250)
  DPoutput = DPmcmc(
    iterations = niter,
    thin = 25,
    data$YList,
    data$TList,
    data$ZList,
    Beta_,
    Lambda2_,
    Nulam_,
    C_,
    Tau2_,
    Nutau_,
    Theta_,
    data$idstart,
    data$idend,
    SigTheta = sigTheta,
    a = 0.1
  )
  time5 = Sys.time()-tick
  
  save(data, DFPoutput, DPoutput, time4, time5, file = paste0("DFPSimOutput/chain",chain_i,".RData"))
  
  return(NULL)
}
stopCluster(cl) # Terminate parallelization