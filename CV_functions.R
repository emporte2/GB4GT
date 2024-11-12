###############################################################################
# R function: this function performs cross-validation for gradient boosting
#             for group testing data using regression trees as the weak learner.
# Input: 
# maxdepth = maximum tree depth, to be selected using CV
# minbucket = minimum observations in any terminal node, to be selected using CV
# cp = complexity parameter, fixed in this example
# epoch = maximum number of boosting iterations to perform
# C = matrix identifying subjects assigned to the kth masterpool, format first column is number of individuals, then followed by indices of individuals 
# B = matrix identifying the test results associated with kth masterpool, first column is number of results, followed by indices of test results 
# Y = matrix whose ith row is of the form (col1=0, col2=number of pools ith individual involved in, col3:col(3+col2-1)=pools ith individual was assigned to)
# Yd = starting values for individual-level outcomes
# X = design matrix beginning with a column of 1's corresponding to an intercept
# Z = a matrix of testing responses
# Se = sensitivity
# Sp = specificity
# k.train = indices of subgroups (out of K) used for training
# k.test = indices of subgroups (out of K) used for test/validation
# ind.train = indices of individuals assigned to training set
# ind.test = indices of individuals assigned to test set

fit.model.tree <- function(maxdepth, minbucket, cp, 
                           epoch, 
                           C, B, Y, Yd, X, Z, Se, Sp, 
                           k.train, k.test, ind.train, ind.test){
  
  CV <- rep(NA, epoch)
  
  test <- cbind(Yd[-ind.train], X[-ind.train,-1])
  test <- as.data.frame(test)
  
  p <- rep(mean(Yd[ind.train]),N)
  p0 <- mean(Yd[ind.train])
  p0te <- mean(Yd[ind.train])
  
  F0 <- log(p0/(1-p0))
  F0te <- log(p0te/(1-p0te))
  
  for(i in 1:epoch){
    
    SR <- grad.GT.train(p,C,B,Y,Z,Se,Sp,k.train)
    
    train <- cbind(SR[ind.train], X[ind.train,-1])
    colnames(train)[1] <- "SR"
    train <- as.data.frame(train)
    fit <- rpart(SR~., data=train, method  = "anova", 
                 control = list(maxdepth = maxdepth, minbucket = minbucket, cp=cp))
    
    htr <-  predict(fit, newdata=train)
    hte <-  predict(fit, newdata=test) 
    
    lam <- optimize(loglik.optim,interval=c(-10,10), F0=F0, h=htr, C=C, B=B, N=N, k.ind=k.train, obs.ind=ind.train)$minimum
    
    F0 <- F0 + lam*htr
    F0te <- F0te + lam*hte
    
    p0 <- exp(F0)/(1+exp(F0))
    p0te <- exp(F0te)/(1+exp(F0te))
    
    p[ind.train] <- p0
    p[-ind.train] <- p0te
    
    CV[i] <- loglik(C,B,N,p,k.test)
  }
  return(CV=CV)
}

###############################################################################
# R function: this function performs cross-validation for gradient boosting
#             for group testing data using kernel smoothing as the weak learner.
# Input: 
# bandwidth = kernel smoothing bandwith, to be chosen by CV
# epoch = maximum number of boosting iterations to perform
# C = matrix identifying subjects assigned to the kth masterpool, format first column is number of individuals, then followed by indices of individuals 
# B = matrix identifying the test results associated with kth masterpool, first column is number of results, followed by indices of test results 
# Y = matrix whose ith row is of the form (col1=0, col2=number of pools ith individual involved in, col3:col(3+col2-1)=pools ith individual was assigned to)
# Yd = starting values for individual-level outcomes
# X = design matrix beginning with a column of 1's corresponding to an intercept
# Z = a matrix of testing responses
# Se = sensitivity
# Sp = specificity
# k.train = indices of subgroups (out of K) used for training
# k.test = indices of subgroups (out of K) used for test/validation
# ind.train = indices of individuals assigned to training set
# ind.test = indices of individuals assigned to test set

fit.model.kernel <- function(bandwidth, 
                             epoch, 
                             C, B, Y, Yd, X, Z, Se, Sp, 
                             k.train, k.test, ind.train, ind.test){
  
  CV <- rep(NA, epoch)
  
  test <- cbind(Yd[-ind.train], X[-ind.train,-1])
  colnames(test) <- c("Z","X")
  test <- as.data.frame(test)
  
  p <- rep(mean(Yd[ind.train]),N)
  p0 <- mean(Yd[ind.train])
  p0te <- mean(Yd[ind.train])
  
  F0 <- log(p0/(1-p0))
  F0te <- log(p0te/(1-p0te))
  
  for(i in 1:epoch){
    
    SR <- grad.GT.train(p,C,B,Y,Z,Se,Sp,k.train)
    
    train <- cbind(SR[ind.train], X[ind.train,-1])
    colnames(train) <- c("SR","X")
    train <- as.data.frame(train)
    htr <- ksmooth(train$X, train$SR, kernel="normal",bandwidth=bandwidth)$y    
    hte <- ksmooth(train$X, train$SR, x.points=X[-ind.train,-1], kernel="normal", bandwidth=bandwidth)$y
    
    # learning rate optimized
    lam <- optimize(loglik.optim,interval=c(-10,10), F0=F0, h=htr, C=C, B=B, N=N, k.ind=k.train, obs.ind=ind.train)$minimum
    
    F0 <- F0 + lam*htr
    F0te <- F0te + lam*hte
    
    p0 <- exp(F0)/(1+exp(F0))
    p0te <- exp(F0te)/(1+exp(F0te))
    
    p[ind.train] <- p0
    p[-ind.train] <- p0te
    
    CV[i] <- loglik(C,B,N,p,k.test)
  }
  return(CV=CV)
}

###############################################################################
# R function: this function performs cross-validation for gradient boosting
#             for group testing data using splines as the weak learner.
# Input: 
# nint = number of interior knots, to be chosen by CV
# epoch = maximum number of boosting iterations to perform
# C = matrix identifying subjects assigned to the kth masterpool, format first column is number of individuals, then followed by indices of individuals 
# B = matrix identifying the test results associated with kth masterpool, first column is number of results, followed by indices of test results 
# Y = matrix whose ith row is of the form (col1=0, col2=number of pools ith individual involved in, col3:col(3+col2-1)=pools ith individual was assigned to)
# Yd = starting values for individual-level outcomes
# X = design matrix beginning with a column of 1's corresponding to an intercept
# Z = a matrix of testing responses
# Se = sensitivity
# Sp = specificity
# k.train = indices of subgroups (out of K) used for training
# k.test = indices of subgroups (out of K) used for test/validation
# ind.train = indices of individuals assigned to training set
# ind.test = indices of individuals assigned to test set

fit.model.spline <- function(nint, 
                             lam, epoch, 
                             C, B, Y, Yd, X, Z, Se, Sp, 
                             k.train, k.test, ind.train, ind.test){
  
  CV <- rep(NA, epoch)
  
  test <- cbind(Yd[-ind.train], X[-ind.train,-1])
  colnames(test) <- c("Z","X")
  test <- as.data.frame(test)
  
  p <- rep(mean(Yd[ind.train]),N)
  p0 <- mean(Yd[ind.train])
  p0te <- mean(Yd[ind.train])
  
  F0 <- log(p0/(1-p0))
  F0te <- log(p0te/(1-p0te))
  
  for(i in 1:epoch){
    
    SR <- grad.GT.train(p,C,B,Y,Z,Se,Sp,k.train)
    
    train <- cbind(SR[ind.train], X[ind.train,-1])
    colnames(train) <- c("SR","X")
    train <- as.data.frame(train)
    
    knots <- quantile(train$X, prob=seq(0,1,length.out=(nint+2))[-c(1,nint+2)])
    Btrain <- bs(train$X, knots=knots)
    Btest <- bs(test$X, knots=knots)
    
    fit <- lm(train$SR ~ Btrain, data=train)
    
    htr <- cbind(1,Btrain)%*%coef(fit)
    hte <- cbind(1,Btest)%*%coef(fit)
    
    F0 <- F0 + lam*htr
    F0te <- F0te + lam*hte
    
    p0 <- exp(F0)/(1+exp(F0))
    p0te <- exp(F0te)/(1+exp(F0te))
    
    p[ind.train] <- p0
    p[-ind.train] <- p0te
    
    CV[i] <- loglik(C,B,N,p,k.test)
  }
  return(CV=CV)
}
