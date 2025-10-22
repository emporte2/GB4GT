###############################################################################
# Read in necessary functions to run the example. 
# Set R directory to your location that conatins the GB4GT R scripts.

source("GB_functions.R")
library(rpart)
library(splines)

# for illustration, we will specify tuning parameter values for the weak 
# learners in this example.  To obtain complete cross-validation results, you may 
# use the functions in CV_functions.R to select your desired tuning parameters.
# source("CV_functions.R")

###############################################################################
# Generate group testing data set according to Dorfman testing protocol.

set.seed(449)

N <- 10000  # sample size
cj <- 4     # pool size
Se <- 0.95  # sensitivity
Sp <- 0.98  # specificity

X <- cbind(rep(1,N),sort(runif(N,0,2*pi))) # design matrix
p.true <- 1/(1+exp(-(-2+0.5*sin(X[,2]))))  # true infection probabilities
Y.true <- rbinom(N,1,p.true)               # true infection statuses

# Simulate group testing (Dorfman testing)
test.res <- Dorfman.decode.same.error(Y.true,Se,Sp,cj)

Z <- test.res$Z # A matrix of testing responses whose jth row is of the form 
#               (col1=Z_j, col2=cj, col3=assay used, col4:col(4+cj-1)=indices 
#               of the individuals assigned to the jth pool) 

Y <- test.res$Y # matrix whose ith row is of the form (col1=0, col2=number of 
#                 pools ith individual involved in, col3:col(3+col2-1)=pools ith
#                 individual was assigned to)

# Build the C matrix for tracking subjects assigned to the kth masterpool
K<-N/cj
C<-matrix(NA, nrow=K, ncol=cj+1)
for(k in 1:K){
  C[k,]<-c(cj, ((k-1)*cj+1):(k*cj))
}

# Build the B matrix for tracking tests associated with kth masterpool

# Note using 500 here to make sure we have enough room, it is trimmed off below
B<-matrix(NA,nrow= K, ncol=500)  

for(k in 1:K){
  
  mps<-C[k,1]
  ind<-C[k,2:(mps+1)]
  temp<-NULL
  
  for(i in ind){
    np<-Y[i,2]
    temp<-c(temp,Y[i,3:(np+2)])
  }
  temp<-unique(temp)
  nt<-length(temp)
  B[k,1]<-nt
  B[k,2:(nt+1)]<-temp
}
B<-B[,apply(is.na(B),2,mean)<1]   # Trimmed here

Yd<-rep(-99,N)
for(i in 1:dim(Y)[1]){
  tests<-Y[i,3:(Y[i,2]+2)]  # identifies the pool(s) that individual i is in
  Yd[i]<-prod(Z[tests,1])   # multiplies observed (group) test outcomes that individual i was involved in
}                           # Yd indicates if an individual ever tested positive

################################################################################
# Perform gradient boosting to obtain estimated infection probabilities for
# all N individuals using *regression trees* as the weak learner

# initialize values for gradient boosting
phat <- mean(Yd)
F0 <- rep(log(phat/(1-phat)),N)
p0 <- exp(F0)/(1+exp(F0))

epoch <- 10
maxdepth <- 1
minbucket <- 300
cp <- 1e-16

for(i in 1:epoch){
  
  SR <- grad.GT(p0,C,B,Y,Z,Se,Sp)
  
  train <- cbind(SR, X[,-1])
  colnames(train) <- c("Z","X")
  train <- as.data.frame(train)
  
  fit <- rpart(Z~X, data=train, method  = "anova", 
               control = list(maxdepth = maxdepth,
                              minbucket = minbucket,
                              cp=cp))
  h <-  predict(fit, newdata=train)
  
  lam <- optimize(loglik.optim,interval=c(-10,10), F0=F0, h=h, C=C, B=B, N=N, k.ind=1:K, obs.ind=1:N)$minimum
  
  F0 <- F0 + lam*h
  
  p0 <- exp(F0)/(1+exp(F0))
}

F0.trees <- F0

################################################################################
# Perform gradient boosting to obtain estimated infection probabilities for
# all N individuals using *kernel smoothing* as the weak learner

# initialize values for gradient boosting
phat <- mean(Yd)
F0 <- rep(log(phat/(1-phat)),N)
p0 <- exp(F0)/(1+exp(F0))

epoch <- 5
bandwidth <- 2

for(i in 1:epoch){
  
  SR <- grad.GT(p0,C,B,Y,Z,Se,Sp)
  
  train <- cbind(SR, X[,-1])
  colnames(train) <- c("SR","X")
  train <- as.data.frame(train)
  
  h <- ksmooth(train$X, train$SR, kernel="normal", bandwidth=bandwidth)$y  

  lam <- optimize(loglik.optim,interval=c(-10,10), F0=F0, h=h, C=C, B=B, N=N, k.ind=1:K, obs.ind=1:N)$minimum
  
  F0 <- F0 + lam*h
  
  p0 <- exp(F0)/(1+exp(F0))
}

F0.ks <- F0

################################################################################
# Perform gradient boosting to obtain estimated infection probabilities for
# all N individuals using *splines* as the weak learner

# initialize values for gradient boosting
phat <- mean(Yd)
F0 <- rep(log(phat/(1-phat)),N)
p0 <- exp(F0)/(1+exp(F0))

epoch <- 2
nint <- 1

for(i in 1:epoch){
  
  SR <- grad.GT(p0,C,B,Y,Z,Se,Sp)
  
  train <- cbind(SR, X[,-1])
  colnames(train) <- c("SR","X")
  train <- as.data.frame(train)
  
  knots <- quantile(train$X, prob=seq(0,1,length.out=(nint+2))[-c(1,nint+2)])
  Btrain <- bs(train$X, knots=knots)
  
  fit <- lm(train$SR ~ Btrain, data=train)
  
  h <- cbind(1,Btrain)%*%coef(fit)
  
  lam <- optimize(loglik.optim,interval=c(-10,10), F0=F0, h=h, C=C, B=B, N=N, k.ind=1:K, obs.ind=1:N)$minimum
  
  F0 <- F0 + lam*h
  
  p0 <- exp(F0)/(1+exp(F0))
}

F0.splines <- F0

###############################################################################
# plot results

par(mfrow=c(1,3))
plot(X[,2], 0.5*sin(X[,2])-2, col="azure3", ylim=c(-2.8,-1.4),
     xlab="x", ylab="F(x)", main="GB with \n regression trees")
lines(X[,2], F0.trees, col="black")

plot(X[,2], 0.5*sin(X[,2])-2, col="azure3", ylim=c(-2.8,-1.4),
     xlab="x", ylab="F(x)", main="GB with \n kernel smoothing")
lines(X[,2], F0.ks, col="black")

plot(X[,2], 0.5*sin(X[,2])-2, col="azure3", ylim=c(-2.8,-1.4),
     xlab="x", ylab="F(x)", main="GB with \n splines")
lines(X[,2], F0.splines, col="black")


