###################################

# Multivariate normal distribution

###################################

library(MASS)
X<-fractions(c(1/2,1,3,1,3,1/4,2,-1,6,5,6,-1/2))
X
A<-fractions(matrix(X,nrow = 4,ncol = 3))
A

############################################################

# Spectral decomposition method for generating Nd(µ, Σ) samples

###############################################################

# mean and covariance parameters
mu <- c(0, 0)
Sigma <- matrix(c(1, .9, .9, 1), nrow = 2, ncol = 2)
Sigma

# The eigen function returns the eigenvalues and eigenvectors of a matrix

rmvn.eigen <-
  function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    ev <- eigen(Sigma, symmetric = TRUE)
    lambda <- ev$values
    V <- ev$vectors
    R <- V %*% diag(sqrt(lambda)) %*% t(V)
    Z <- matrix(rnorm(n*d), nrow = n, ncol = d)
    X <- Z %*% R + matrix(mu, n, d, byrow = TRUE)
    X
  }

# generate the sample
X <- rmvn.eigen(1000, mu, Sigma)
plot(X, xlab = "x", ylab = "y", pch = 20)
print(colMeans(X))
print(cor(X))

###############################################################
# SVD Method of generating Nd(µ, Σ) samples
#######################################################

rmvn.svd <-
  function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    S <- svd(Sigma)
    R <- S$u %*% diag(sqrt(S$d)) %*% t(S$v) #sq. root Sigma
    Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
    X <- Z %*% R + matrix(mu, n, d, byrow=TRUE)
    X
  }

##########################################################
# Choleski factorization method of generating Nd(µ, Σ) samples
###########################################################

rmvn.Choleski <-
  function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    Q <- chol(Sigma) # Choleski factorization of Sigma
    Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
    X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
    X
  }
y <- subset(x=iris, Species=="virginica")[, 1:4]
mu <- colMeans(y)
Sigma <- cov(y)
mu
Sigma
#now generate MVN data with this mean and covariance
X <- rmvn.Choleski(200, mu, Sigma)
pairs(X)
View(iris)

##############################################################
# Comparing Performance of Generators
###############################################################

library(MASS)
library(mvtnorm)
n <- 100 #sample size
d <- 30 #dimension
N <- 2000 #iterations
mu <- numeric(d)
set.seed(100)
system.time(for (i in 1:N)
  rmvn.eigen(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
  rmvn.svd(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
  rmvn.Choleski(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
  mvrnorm(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
  rmvnorm(n, mu, cov(matrix(rnorm(n*d), n, d))))
set.seed(100)
system.time(for (i in 1:N)
  cov(matrix(rnorm(n*d), n, d)))

###########################################################
# Monte Carlo integration
##########################################################

m<-10^6
X<-runif(m)
2*exp(-6)*mean(exp(-4*X)) 
###############################################
m<-10^6
y<-runif(m)
y
0.5-2.5/sqrt(2*pi)*mean(exp(-1/2*(2.5*y)^2)) 
######################################################
# Find Q(-2.5)
m<-10^5
X<-runif(m,-2.5,0)
exp(-1/2*X^2)
0.5-2.5/sqrt(2*pi)*mean(exp(-1/2*X^2))
