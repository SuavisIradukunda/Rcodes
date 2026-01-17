# FOR STUDENT defining degree of freedom for t-distribution
n_df<-10
random_t<-rt(15, df = n_df)
print(random_t)

# FOR FISHER defining degree of freedom for F-distribution
n_df_F <- 10
m_df_F<-15
random_f<-rf(15,df1 = n_df_F, df2 = m_df_F)
print(random_f)

##########################################################

# Generate Student t using n-fold convolution

##########################################################
rt_convolution <- function(n, size = 1000) {
  t_vals <- numeric(size)
  
  for (i in 1:size) {
    # Step 1: Generate n standard normals
    Z  <- rnorm(1)
    Zi <- rnorm(n)
    V  <- sum(Zi^2) # n-fold convolution 
    t_vals[i] <- Z / sqrt(V / n)
  }
  return(t_vals)
}

# Simulation
set.seed(123)
x <- rt_convolution(n = 10, size = 1000)

# Histogram and theoretical density
hist(x, probability = TRUE,
     main = "Student t(10) via n-fold Convolution",
     xlab = "x", col = "lightgray", border = "white")

curve(dt(x, df = 10),
      add = TRUE, col = "red", lwd = 2)



#####################################################

# Generate Fisher F using n-fold convolution

###################################################
rf_convolution <- function(n, m, size = 1000) {
  f_vals <- numeric(size)
  
  for (i in 1:size) {
    U <- sum(rnorm(n)^2)
    V <- sum(rnorm(m)^2)
    f_vals[i] <- (U / n) / (V / m)
  }
  return(f_vals)
}

# Simulation
set.seed(123)
y <- rf_convolution(n = 5, m = 10, size = 1000)

# Histogram and theoretical density
hist(y, probability = TRUE,
     main = "F(5,10) via n-fold Convolution",
     xlab = "y", xlim = c(0, 5),
     col = "lightgray", border = "white")

curve(df(x, df1 = 5, df2 = 10),
      add = TRUE, col = "blue", lwd = 2)

