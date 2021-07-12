
library(MASS)

# Setup
n <- 2000
L <- 21 # number of candidate instruments
s <- 9  # number of invalid instruments

gamma1 <- runif(L,1.5,2.5);
gamma2 <- runif(L,1.5,2.5);
gamma <- cbind(gamma1, gamma2) # first stage coefficients
cgamma <- 0.4

# direct effects of instruments on the outcome
pi <- c(rep(1,s),rep(0,L-s)) 
cpi <- 0.6

beta <- c(0.5, 1) # causal effects

# variance-covariance matrix of instruments
Sigmaz <- matrix(rep(0,L*L), nrow=L)
for(i in 1:L){
  for(j in 1:L){
    Sigmaz[i,j] <- 0.5^abs(i-j) 
  }
}

# variance and covariance of the error terms
epsvar <- 1
rho1 <- 0.25
rho2 <- 0.3

# Generate data
error <- rnorm(n,0,epsvar)
Z <- matrix(mvrnorm(n,mu=rep(0,L),Sigma=Sigmaz), ncol=L)

errorD <- matrix(rnorm(2*n,0,1),nrow=n, ncol=2) + cbind(rho1*error, rho2*error)
D <- 0.1 + Z %*% gamma * cgamma + errorD
colnames(D) = NULL

# two exogenous control variables
X <-  matrix(mvrnorm(n,mu=rep(0,2),Sigma=diag(1,2)), ncol=2) 
eta <- c(0.5, 0.8)

Y <- 0.1  + D %*% beta + X %*% eta + Z %*% pi * cpi + error

# IV selection and estimation.
source("MVMR Adptive Lasso Functions.R")
output.cv <- MVadap.cv(Y,D,Z,X)
output.dt <- MVadap.dt(Y,D,Z,X)
