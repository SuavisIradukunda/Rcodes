
# Read Notes on Christian Robert, p.67-71, and 
# Maria L. Rizzo, p.69-71
# Optimization of beta(2.7,6.3)/beta(2,6) on (0,1)
  optimize(f=function(x){dbeta(x,2.7,6.3)/dbeta(x,2,6)},
          maximum=T,interval=c(0,1))$objective
#
  Nsim=2500
  a=2.7;b=6.3
  M=1.671808

#uniform over (0,M)
u=runif(Nsim,max=M)

#generation from g=beta(2,6), the instrumental density
y=rbeta(Nsim,2,6)

# Generate a random sample beta(2.7,6.3) of the targeted 
# from instrumental beta(2,6) by accept-reject method
  x=y[u<dbeta(y,a,b)]
  hist(x,breaks = 30  ,probability = T)
  y<-seq(0,1,0.01)
  lines(y,dbeta(y,2.7,6.3),col="red")
  
  
####################################################
 # EXERCISE/ASSIGNMENT FOR TOMMORROW: Christian Robert,
 # Exercise 2.8., page 72
####################################################
   
  # Mixtures
  
  n <- 1000
  nu <- 2
  X <- matrix(rnorm(n*nu), n, nu)^2 #matrix of sq. normals
  #sum the squared normals across each row: method 1
  y <- rowSums(X)
  hist(y,probability = T,100)
  X<-seq(0,20,0.01)
  lines()
  #method 2
  y <- apply(X, MARGIN=1, FUN=sum) #a vector length n
  
  ###############################################
  
  NS<-1000
  X1<-rnorm(NS,0,1)
  X2<-rnorm(NS,3,1)
  X<X1+X2
  hist(X,probability = T,200)
  mean(X)
  sqrt(var(X))
