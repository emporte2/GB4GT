###############################################################################
######### Functions for performing gradient boosting on pooled data ###########

# R function: this function computes subject-specific gradients for all observations
# Input: 
# p = length N vector of probabilities
# C = matrix identifying subjects assigned to the kth masterpool, format first column is number of individuals, then followed by indices of individuals 
# B = matrix identifying the test results associated with kth masterpool, first column is number of results, followed by indices of test results 
# Y = matrix whose ith row is of the form (col1=0, col2=number of pools ith individual involved in, col3:col(3+col2-1)=pools ith individual was assigned to)
# Z = a matrix of testing responses
# Se = sensitivity
# Sp = specificity

grad.GT<-function(p,C,B,Y,Z,Se,Sp){
  
  N<-dim(Y)[1]  # number of (individual) observations
  K<-dim(C)[1]  # number of non-overlapping subgroups for forming the likelihood
  
  dprob.Z<-rep(-99,N)
  for(k in 1:K){
    
    ind<-C[k,2:(C[k,1]+1)]  #Pull individuals associated with kth group
    ptemp<-p[ind]           #Pull probs for individuals associated with kth group
    cj<-length(ind)         # number of individuals in group k
    tmat<-expand.grid(replicate(cj, 0:1, simplify = FALSE))  # all 2^cj possible combinations of 0's and 1's 
    prod.pi<-apply(t(t(tmat)*ptemp) + t(t(abs(1-tmat))*(1-ptemp)),1,prod)  # gets a probability for each combinations of 0's and 1's for th kth group
    
    SeSp<-1
    # Zd for test outcome, Zt for true outcome
    for(t in 1:B[k,1]){
      Zd<-Z[B[k,(t+1)],1]  # test outcomes for kth group
      Yid<-Z[B[k,(t+1)],4:(Z[B[k,(t+1)],2]+3)]  # individuals contributing to test outcome
      id<-match(Yid, ind)
      Zt<-apply(as.matrix(tmat[,id]),1,sum)>0  # vector of all true statuses of Zd 
      SeSp<- (Zt*Zd*Se + Zt*(1-Zd)*(1-Se) + (1-Zt)*Zd*(1-Sp) +(1-Zt)*(1-Zd)*Sp)*SeSp
    }
    
    prob.Z<-sum(SeSp*prod.pi)
    
    for(i in ind){
      id<-match(i, ind)
      term<-(1-p[i])*tmat[,id] -p[i]*(1-tmat[,id])
      dprob.Z[i]<-sum(SeSp*prod.pi*term)/prob.Z
    }
  }
  return(dprob.Z)
}

###############################################################################
# R function: this function is used to select the
#             learning rate at each epoch/iteration by optimizing the
#             log-likelihood given the current model fit.
# Input: 
# lam = learning rate
# F0 = current model fit
# h = current weak learner fit
# C = matrix identifying subjects assigned to the kth masterpool, format first column is number of individuals, then followed by indices of individuals 
# B = matrix identifying the test results associated with kth masterpool, first column is number of results, followed by indices of test results 
# N = total number of individuals
# k.ind = indices of subgroups used to calculate likelihood (set to 1:K for entire data set)
# obs.ind = indices of individuals used to calculate likelihood (set to 1:N for entire data set)

loglik.optim <- function(lam,F0,h,C,B,N,k.ind,obs.ind){
  
  F <- F0+lam*h
  
  p.sub <- exp(F)/(1+exp(F))
  p <- rep(NA,N)
  p[obs.ind] <- p.sub
  
  prob.Z <- rep(NA,length(k.ind))
  for(k in k.ind){
    
    ind <- C[k,2:(C[k,1]+1)]  #Pull individuals associated with kth group
    ptemp <- p[ind]           #Pull probs for individuals associated with kth group
    cj <- length(ind)         # number of individuals in group k
    tmat <- expand.grid(replicate(cj, 0:1, simplify = FALSE))  # all 2^cj possible combinations of 0's and 1's 
    prod.pi <- apply(t(t(tmat)*ptemp) + t(t(abs(1-tmat))*(1-ptemp)),1,prod)  # gets a probability for each combinations of 0's and 1's for th kth group
    
    SeSp <- 1
    # Zd for test outcome, Zt for true outcome
    for(t in 1:B[k,1]){
      Zd <- Z[B[k,(t+1)],1]  # test outcomes for kth group
      Yid <- Z[B[k,(t+1)],4:(Z[B[k,(t+1)],2]+3)]  # individuals contributing to test outcome
      id <- match(Yid, ind)
      Zt <- apply(as.matrix(tmat[,id]),1,sum)>0  # vector of all true statuses of Zd 
      SeSp <- (Zt*Zd*Se + Zt*(1-Zd)*(1-Se) + (1-Zt)*Zd*(1-Sp) +(1-Zt)*(1-Zd)*Sp)*SeSp
    }
    
    prob.Z[k] <- sum(SeSp*prod.pi)
  }
  
  loglik <- -sum(log(na.omit(prob.Z)))
  
  return(loglik=loglik)
}

###############################################################################
# R function: This function simulates Dorfman decoding and stores the 
#             testing responses in accordance to the data structure 
#             required to perform gradient boosting for group testing
#             data as described in Porter et al. (202X)
#
# Input: 
# Y.true = The true statuses of the individuals
# Se = The testing sensitivity used for both pools and individual testing
# Sp = The testing specificity used for both pools and individual testing
# cj = pool size to be used

Dorfman.decode.same.error<-function(Y.true,Se,Sp,cj){
  N<-length(Y.true)
  Jmax<-N+N/cj 
  J<-1
  
  Y<-matrix(-99,nrow=N, ncol=4) 
  Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 
  
  
  for(j in 1:(N/cj)){
    prob<-ifelse(sum(Y.true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
    Z[J,1]<-rbinom(n=1,size=1,prob=prob)
    Z[J,2]<-cj
    Z[J,3]<-1
    Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
    Y[((j-1)*cj+1):(j*cj),2]<-1
    Y[((j-1)*cj+1):(j*cj),3]<-J
    J<-J+1
    if(Z[J-1,1]==1){
      for(k in ((j-1)*cj+1):(j*cj)){
        prob<-ifelse(Y.true[k]>0,Se,1-Sp)
        Z[J,1]<- rbinom(n=1,size=1,prob=prob)
        Z[J,2]<-1
        Z[J,3]<-1
        Z[J,4]<-k
        Y[k,2]<-2
        Y[k,4]<-J
        J<-J+1
      }
    }
  }
  
  J<-J-1
  Z<-Z[1:J,]
  
  return(list("Z"=Z, "Y"=Y))
}

