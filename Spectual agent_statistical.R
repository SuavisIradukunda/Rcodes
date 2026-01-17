
# Set seed for reproducibility
set.seed(123)

# Mean vector
mu <- c(0, 1, 2)

# Covariance matrix
Sigma <- matrix(c( 1.0, -0.5,  0.5,
                   -0.5,  1.0, -0.5,
                   0.5, -0.5,  1.0),
                nrow = 3, byrow = TRUE)

# Step 1: Spectral (eigen) decomposition
eig <- eigen(Sigma)

Q <- eig$vectors                 # Eigenvectors
Lambda_sqrt <- diag(sqrt(eig$values))  # Square roots of eigenvalues

# Step 2: Generate standard normal variables
n <- 200
Z <- matrix(rnorm(n * 3), nrow = 3)

# Step 3: Construct multivariate normal sample
X <- Q %*% Lambda_sqrt %*% Z
X <- t(X) + mu   # transpose and add mean

# Convert to data frame
X <- as.data.frame(X)
colnames(X) <- c("X1", "X2", "X3")

# View first rows
head(X)

