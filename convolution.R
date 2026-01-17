
rt_convolution <- function(n, size = 1) {
  t_values <- numeric(size)
  
  for (k in 1:size) {
    # Step 1: Generate Z ~ N(0,1)
    Z <- rnorm(1)
    
    # Step 2: Generate n independent N(0,1) variables
    Zi <- rnorm(n)
    
    # Step 3: n-fold convolution (chi-square)
    V <- sum(Zi^2)
    
    # Step 4: Student's t variable
    t_values[k] <- Z / sqrt(V / n)
  }
  
  return(t_values)
}
rt(1000, df = 10)   # Built-in Student t

set.seed(123)
x <- rt_convolution(n = 10, size = 1000)

hist(x, probability = TRUE, main = "Student t(10) via Convolution")
curve(dt(x, df = 10), add = TRUE, col = "red", lwd = 2)

rf_convolution <- function(n, m, size = 1) {
  f_values <- numeric(size)
  
  for (k in 1:size) {
    # Step 1: Generate n standard normals
    Xi <- rnorm(n)
    U <- sum(Xi^2)   # n-fold convolution
    
    # Step 2: Generate m standard normals
    Yj <- rnorm(m)
    V <- sum(Yj^2)   # m-fold convolution
    
    # Step 3: Fisher F variable
    f_values[k] <- (U / n) / (V / m)
  }
  
  return(f_values)
}

set.seed(123)
y <- rf_convolution(n = 5, m = 10, size = 1000)

hist(y, probability = TRUE, main = "F(5,10) via Convolution", xlim = c(0, 5))
curve(df(x, df1 = 5, df2 = 10), add = TRUE, col = "blue", lwd = 2)


rf(1000, df1 = 5, df2 = 10)  # Built-in F

