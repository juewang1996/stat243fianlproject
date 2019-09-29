#Test if the x in the updated matrix is in an ascending order
update_mat <- function(mat,x_star,fun,deriv,a,b){
  index <- sum(mat[,1] < x_star) #find the location to insert x*
  tmp <- matrix(NA,1,3)
  mat <- rbind(mat,tmp)
  mat[(index+2):nrow(mat),] <- mat[(index+1):(nrow(mat)-1),]
  mat[(index+1),] <- c(x_star,fun(x_star),deriv(x_star,fun,a,b))
  return(mat)
}
#test using the normal distribution
h <- function(x){
  return(log(g(x)))
}
g <- function(x){
  return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
}
Deriv_sub <- function(x, fun=h, a, b){
  if (x == a) {return ((fun(x + 1e-9)-fun(x))/1e-9)}
  if (x == b) {return ((fun(x) - fun(x - 1e-9))/1e-9)}
  if (a <= x && x <= b) {return((fun(x + 1e-9)-fun(x - 1e-9))/2e-9)}
}
mat_test <- matrix(c(-2, 0, 2, h(-2),h(0),h(2),Deriv_sub(-2,h,-2,2),Deriv_sub(0,h,-2,2),Deriv_sub(2,h,-2,2)),
              ncol = 3)
mat_update <- update_mat(mat_test,1,h,Deriv_sub,-2,2)
expect_equal(mat_update[,1], sort(c(1,mat_test[,1])))
if(sum(expect_equal(mat_update[,1], sort(c(1,mat_test[,1]))))){ 
  print("the x in the updated matrix is in an ascending orderl, update_mat unit test passed")
  }

#Test if the sample is from the density we interested
sample_x <- function(z, mat, n=1){
  uk1 <- exp(mat[,2]+(z[-1]-mat[,1])*mat[,3])
  uk2 <- exp(mat[,2]+(z[-length(z)]-mat[,1])*mat[,3])
  p <- rep(NA,nrow(mat))
  p <- (uk1-uk2)/mat[,3]
  
  # normalize p
  p_norm <- p/sum(p)
  w <- runif(n)
  # determine the range we want to sample from
  i <- sum(cumsum(p_norm) < w) + 1
  
  # sample x using inverse cdf
  samp_x <- (log(p[i]*mat[i,3]*runif(n)+exp(mat[i,2]+(z[i]-mat[i,1])*mat[i,3]))-mat[i,2])/mat[i,3]+mat[i,1]
  return(samp_x)
}

set.seed(5) 
a <- 0 
b <- 1 
lambda <- 1
mat_test_2 <- matrix(c(a, b,log(dexp(a, rate=lambda)), log(dexp(b, rate= lambda)),-lambda,-lambda), nrow=2)
z <- c(0, 0.5, Inf)
sample <- replicate(1000,sample_x(z,mat_test_2,1))
out.test <- ks.test(sample, pexp)
if(sum(expect_gt(out.test$p.value, 0.05))){ 
  print("sample matches exponential, sample_x unit test passed") } else{
  print("sample does not match exponential, sample_ex unit test failed") 
  }
expect_gt(out.test$p.value, 0.05)
