gendata2 = function(mEta = 0, myseed = 250, N = 50, J = 5, dT = 5, dZ = 3, dX = 3, a = 0.5, ad = 0.5){
  
  # mEta = -3
  # myseed = 36
  # N = 50 # number of subjects
  # J = 5 # number of periods
  # dT = 5 # col dimension of spline basis functions
  # dZ = 3 # col dimension of global effect matrices
  # dX = 3 # dimension of reallocation coefficients
  # a=ad=0.5
  
  M = matrix(30,nrow=N,ncol=J) # number of timepoints per subject within each period
  
  # --- todo: change this such that each i share the same global covariate
  set.seed(myseed)
  ZList = list() # each element is a sum_it x dZ matrix
  for( j in 1:J ){
    Zj = NULL
    for( i in 1:N ){
      Zij0 = rnorm(dZ)
      Zij = matrix(rep(Zij0, M[i,j]), ncol = dZ, byrow = TRUE)
      Zj = rbind(Zj, Zij)
    }
    ZList[[j]] = Zj
  }
  Theta = matrix(rep(c(0.5, -0.5, 0.3),J), ncol=J) 
  # Theta = matrix(rep(c(0, 0, 0),J), ncol=J) 
  
  # simulate obs timepoint obsTime and make spline basis TList
  # each element of obsTime size: sum_it x J
  library(splines)
  
  TList = list()
  obsTime = list()
  for( j in 1:J ){
    ot_j = NULL # obs time point for all i at period j
    for( i in 1:N ){
      ot = sort(runif(M[i,j])) # obs time point for i at period j
      ot_j = c(ot_j,ot)
    }
    obsTime[[j]] = ot_j # using a list since the number of observation for different j might not be the same
    T_j = bs(obsTime[[j]], df = dT, intercept = TRUE)
    TList[[j]] = T_j
  }
  
  # set.seed(125)
  # xx = sort(runif(50))
  # x = bs(xx, df = 5, intercept = T)
  # plot(xx, x[,1], type='l')
  # lines(xx, x[,2], col = 2)
  # lines(xx, x[,3], col = 3)
  # lines(xx, x[,4], col = 4)
  # lines(xx, x[,5], col = 5)
  
  # Simulate some covariates to affect the clusters transitioning
  
  Eta = rnorm(dX, mEta, 0.5)
  # Eta = rep(2,dX)
  
  X = array( 0, dim = c(N, dX, J) )
  for(j in 1:J){
    xtmp = matrix(rnorm(N*(dX), mean = 0, sd = 0.5), nrow = N)
    # xtmp = matrix(0, nrow = N, ncol = dX)
    xtmp[,1] = 1
    X[,,j] = xtmp
  }
  
  # strength of temporal dependence
  Amat = matrix(NA, N, J) 
  for(j in 1:J){
    Amat[,j] = (X[,,j]) %*% Eta
  }
  # rowMeans(Amat) # matrix of X*Eta
  # sum(rowMeans(Amat))
  alpMat = exp(Amat) / ( 1 + exp(Amat) )
  colMeans(alpMat)
  
  # generate the moving indicator
  
  gMat <- matrix(1, nrow=J, ncol=N)
  gMat[1,]=0
  for(j in 2:J){
    gMat[j,] = rbinom(N,1,alpMat[,j])
  }
  
  # get the start/end position of observations
  IDstart = apply(M, MARGIN = 2, cumsum)
  IDstart = rbind(0, IDstart)
  IDstart = IDstart[-nrow(IDstart),]
  IDend = apply(M, MARGIN = 2, cumsum)-1
  
  # From Page's paper
  # This function was used to generate data when exploring correlation on the 
  # data level.
  # March 17 2018
  # N - number of experimental units
  # a - DP scale (concentration) parameter
  # alpha - temporal dependence parameter (proportion of being fixed)
  # ntime - number of time points
  # hsMean - mean associated with atom generation
  # hsSig - sd associated with atom generation
  # sig - sd associated with data generation (NOT USED FOR NOW)
  # FirstPart - Allows me to provide the first partition if so desired
  # Type - type of model for the cluster-specific paramters
  #		Random - cluster-specific parameters are drawn randomly across time
  #		Deterministic - cluster-specific parameters are same over time
  #		AR1 - cluster-specific parameters are drawn from an AR(1) process
  #				If this is selected, phi0 and phi1 must be supplied.
  
  ## possible varying effects
  x = seq(0,1,length.out=30)
  f11 = function(x){4*sin(3*x)-2}
  f12 = function(x){-3*sin(3*x)+1.5}
  f13 = function(x){3*cos(3*x)-0.5}
  # f14 = function(x){-1.2*x-1.8}
  f14 = function(x){-3*cos(3*x)+0.5}
  f15 = function(x){3*x}
  f16 = function(x){3*(x-1)^2-1}
  plot(x,f11(x),type='l',ylim = c(-5,5))
  lines(x,f12(x),type='l', col = 2)
  lines(x,f13(x),type='l', col = 3)
  lines(x,f14(x),type='l', col = 4)
  lines(x,f15(x),type='l', col = 5)
  lines(x,f16(x),type='l', col = 6)
  
  f1 = c(f11,f12,f13,f14,f15,f16) 
  
  N <- N;
  # a <- 0.5;
  # ad = 0.5 
  ntime <- J;
  FirstPart = NULL;
  
  library(MASS)
  
  ci <- FirstPart
  if(is.null(FirstPart)){
    ci = rep(1, N)
    ci[sample( 1:N, floor(N/2) )] = 2 ## split subjects to 2 clusters
    K = length(unique(ci))
    Tm = ntime
  }
  
  di = ci # dish
  
  f1tmp = NULL
  for(ddi in di){
    f1tmp = c(f1tmp, f1[[ddi]])
  }

  ciMat <- matrix(NA, nrow=Tm, ncol=N)
  ciMat[1,]  <- ci # for saving the true clustering assignment
  
  diMat <- matrix(0, nrow=Tm, ncol=N)
  diMat[1,]  <- di # for saving the true clustering assignment
  
  PList = list()
  YList = list()
  Psi = numeric(nrow(TList[[1]]))
  Y = numeric(nrow(TList[[1]]))
  
  for(i in 1:N){
    span = IDstart[i,1]:IDend[i,1]
    span = span + 1
    Psi[ span ] = f1tmp[[i]](obsTime[[1]][ span ]) # simulate the outcome
  }
  Psi = Psi + ZList[[1]] %*% Theta[,1]
  Y = rbinom(length(Psi), 1, exp( Psi )/( 1 + exp( Psi ) ) )
  PList[[1]] = Psi
  YList[[1]] = Y 
  
  for(t in 2:Tm){
    r = gMat[t,]
    dnk = which(r == 0) # which subjects are free to move
    fix_subj = which(r == 1) # which subjects are fixed
    
    if(length(dnk) > 0){
      ci[dnk]<- 0
    }	
    
    mh <- tabulate(ci[ci!=0]); # count the existing clusters
    if( (K - length(mh)) > 0 ){ # ?
      mh <- c(mh, rep(0, K-length(mh)))
    }
    
    K <- length(unique(ci[ci!=0])) # number of unique clusters with fixed units
    
    if(sum(r) > 0){ ## reset cluster indicator for the fixed subjects to make sure they start from 1
      ci[ci!=0] <- as.numeric( 
        factor(  ci[ci!=0], labels=1:length( unique(ci[ci!=0]) ) )     
      )
    }
    
    for(k in dnk){
      p <- 1
      if(K>0) p <- c(mh[mh!=0]/(sum(mh[mh!=0])+a), a/(sum(mh[mh!=0])+a))
      
      ci[k] <- sample(1:(K+1), 1, prob=p)
      mh <- table(ci[ci!=0])
      K <- length(unique(ci[ci!=0]))
    } # sample new clusters iteratively for free subjects
    
    ciMat[t,] <- as.numeric(
      factor(
        ci, labels=1:length( unique(ci) )
      )
    )
    
    di = ciMat[t,]
    attn_pos = which(diMat[t-1,fix_subj] != di[fix_subj])
    for(attn_sub in fix_subj[attn_pos]){
      same_tab = which(ciMat[t,] == ciMat[t,attn_sub])
      replace_label = which(di == diMat[t-1,attn_sub])
      old_label = di[attn_sub]
      di[same_tab] = diMat[t-1,attn_sub]
      if(length(replace_label)>0){
        di[replace_label] = old_label
      }
    }

    attn_pos = which(di > 6)
    di[attn_pos] = 1
    gMat[t,attn_pos] = 0
    
    diMat[t,] = di
    
    fjtmp = NULL
    for(ddi in diMat[t,]){
      fjtmp = c(fjtmp, f1[[ddi]])
    }
    
    
    Psi = numeric(nrow(TList[[t]]))
    Y = numeric(nrow(TList[[t]]))
    
    for(i in 1:N){
      span = IDstart[i,t]:IDend[i,t]
      span = span + 1
      Psi[ span ] = fjtmp[[i]]( obsTime[[t]][ span ] )
    }
    Psi = Psi + ZList[[t]] %*% Theta[,t]
    Y = rbinom(length(Psi), 1, exp( Psi )/( 1 + exp( Psi ) ) )
    
    PList[[t]] = Psi
    YList[[t]] = Y
  }
  
  return(list(YList = YList, TList = TList, ZList = ZList, PList = PList, X = X, time = obsTime, Eta = Eta, C = t(ciMat), D = t(diMat), Gamma = t(gMat), idstart = IDstart, idend = IDend))
  
}
