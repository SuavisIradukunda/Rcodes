
# Set parameters
n <- 200
mu <- c(0, 1, 2)

# Define the Matrix (Symmetric)
Sigma <- matrix(c( 1,   -0.5,  0.5, 
                   -0.5,  1,   -0.5, 
                   0.5, -0.5,  1), nrow = 3, ncol = 3)

# Choleski Factorization
L <- t(chol(Sigma))

# Generate 200 observations
set.seed(123)
Z <- matrix(rnorm(3 * n), nrow = 3, ncol = n)
X <- L %*% Z + mu
df <- as.data.frame(t(X))
colnames(df) <- c("X1", "X2", "X3")

# Visualization
pairs(df, 
      main = "Mixed Correlation MVN\n(X1/X3 Positive, others Negative)",
      pch = 19, 
      col = rgb(0.2, 0.6, 0.2, 0.6))



#Spectral decomposition 
# Setup
n <- 200
mu <- c(0, 1, 2)
Sigma <- matrix(c(1, -0.5, 0.5, -0.5, 1, -0.5, 0.5, -0.5, 1), 3, 3)

# Spectral Decomposition
ev <- eigen(Sigma, symmetric = TRUE)
V <- ev$vectors
Lambda_half <- diag(sqrt(ev$values))

# Create transformation matrix A
#Creating the Symmetric Square Root
A_eigen <- V %*% Lambda_half %*% t(V)

# Generate Samples
set.seed(123)
Z <- matrix(rnorm(3 * n), nrow = 3, ncol = n)
X_eigen <- A_eigen %*% Z + mu
df_eigen <- as.data.frame(t(X_eigen))
colnames(df_eigen) <- c("X1", "X2", "X3")

# Plot
pairs(df_eigen, main = "MVN: Spectral Decomposition", col = "steelblue", pch = 19)


# Singular Value Decomposition
s <- svd(Sigma)
U <- s$u
D_half <- diag(sqrt(s$d))

# Create transformation matrix A
# Note: We can use U %*% D_half or U %*% D_half %*% t(V)
A_svd <- U %*% D_half

# Generate Samples
set.seed(123) # Same seed to compare results
Z <- matrix(rnorm(3 * n), nrow = 3, ncol = n)
X_svd <- A_svd %*% Z + mu
df_svd <- as.data.frame(t(X_svd))
colnames(df_svd) <- c("X1", "X2", "X3")

# Plot
pairs(df_svd, main = "MVN: SVD Method", col = "darkblue", pch = 19)
 
