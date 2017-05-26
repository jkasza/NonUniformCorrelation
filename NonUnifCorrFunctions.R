###################################################################
# R code to accompany "Impact of non-uniform correlation structure
#    on sample size and power in multiple-period cluster randomised
#    trials", by J Kasza et al.
# 
# This file contains functions to allow calculation of variances 
# and construction of design matrices 
###################################################################

ExpDecayVar <- function(r,Xmat, m, rho0) {
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  
  T <- ncol(Xmat)
  K <- nrow(Xmat)
  Xvec <- as.vector(t(Xmat))
  
  #Variance matrix for one cluster, with decay in correlation over time
  Vi <- diag(sig2,T) + sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(Xmat)%*%solve(Vi)%*%(matrix(colSums(Xmat),nrow=T, ncol=1))/K )
  return(vartheta)
  
  
}


HooperVar <- function(r, Xmat, m, rho0) {
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  
  T <- ncol(Xmat)
  K <- nrow(Xmat)
  Xvec <- as.vector(t(Xmat))
  
  Vi <-diag(sig2 +(1-r)*sig2CP, T) + matrix(data=sig2CP*r, nrow=T, ncol=T)
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(Xmat)%*%solve(Vi)%*%(matrix(colSums(Xmat),nrow=T, ncol=1))/K )
  return(vartheta)
}


###################################################################
# Construct design matrices
###################################################################

SWdesmat2 <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}


plleldesmat2 <- function(T) {
  if((T-1)%%2 == 0) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T -1)
    Xpllel[1:(T-1)/2,] <- 1
    return(Xpllel)
  }
  
  
  if((T-1)%%2 == 1) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T)
    Xpllel[1:(T)/2,] <- 1
    return(Xpllel)
  }
  
}

pllelbasedesmat2 <- function(T) {
  if((T-1)%%2 == 0) {
    
    Xpllelbase <- matrix(data=0, ncol = T, nrow = T-1)
    Xpllelbase[1:(T-1)/2,2:T] <- 1
    return(Xpllelbase)
  }
  
  
  if((T-1)%%2 == 1) {
    
    Xpllelbase <- matrix(data=0, ncol = T, nrow = T)
    Xpllelbase[1:(T)/2,2:T] <- 1
    return(Xpllelbase)
    
  }
}

crxodesmat2 <- function(T) {
  if((T-1)%%2 == 0) {
    
    Xcrxo <- matrix(data=0, ncol = T, nrow = T-1)
    Xcrxo[1:(T-1)/2, seq(1,T,2)] <- 1
    Xcrxo[((T-1)/2 + 1):(T-1), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  
  
  if((T-1)%%2 == 1) {
    
    Xcrxo <- matrix(data=0, ncol = T, nrow = T)
    Xcrxo[1:(T)/2, seq(1,T,2)] <- 1
    Xcrxo[((T)/2 + 1):(T), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  
}

###################################################################


