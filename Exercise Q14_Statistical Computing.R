library(MASS)

############################################################

# Spectral decomposition method for generating Nd(µ, Σ) samples

###############################################################

# mean and covariance parameters
mu <- c(0,1,2)
Sigma <- matrix(c(1.0, -0.5,  0.5,
                -0.5,  1.0, -0.5,
                0.5, -0.5,  1.0),
nrow = 3, byrow = TRUE)
Sigma
mu

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
X <- rmvn.eigen(200, mu, Sigma)
plot(X, xlab = "x", ylab = "y")
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

# generate the sample
X <- rmvn.eigen(200, mu, Sigma)
plot(X, xlab = "x", ylab = "y")
print(colMeans(X))
print(cor(X))