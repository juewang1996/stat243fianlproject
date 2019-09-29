library(testthat)
library(ars)

#pdf of normal
dnor <- function(x){
  return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
}

#pdf of gamma
dgam <- function(x, theta=2, k=2){
  return((1/(gamma(k)*theta^k))*(x^(k-1))*exp(-x/theta))
}

#pdf of beta(alpha=2,b=2)
dbet <- function(x, alpha=2, b=2){
  return((x^(alpha-1)*(1-x)^(b-1))/beta(alpha, b))
}

#pdf of chi-square distribution
dchi <- function(x) dchisq(x,3)

#pdf of x^3
dx2 <- function(x){
  return(exp(x^2))
}


test_that("empirical distribution matches",{
  #test normal distribution
  print("testing the normal distribution")
  samp_norm <- ars(10000, dnor, a = -10, b = 10)$sample
  ks_norm <- ks.test(samp_norm,rnorm(10000))
  if (expect_that(ks_norm$p.value, function(x){x > 0.05}) == TRUE){
    print("passed with normal distribution")
  } else {print("failed to pass with normal distribution")}

  #test gamma distribution
  print("testing gamma(shape=2,scale=2) distribution")
  samp_gamma <- ars(10000,dgam,a = 0.001, b = 1000)$sample
  ks_gamma <- ks.test(rgamma(10000,shape=2,scale=2),samp_gamma)
  if (expect_that(ks_gamma$p.value, function(x){x > 0.05}) == TRUE){
    print("passed with gamma distribution")
  } else {print("failed to pass with gamma distribution")}

  #test beta distribution
  print("testing beta(2,2) distribution")
  samp_beta <- ars(10000,dbet,a=0.01,b=0.99)$sample
  ks_beta <- ks.test(samp_beta,rbeta(10000,2,2))
  if (expect_that(ks_beta$p.value,function(x){x>0.05}) == TRUE){
    print("passed with beta distribution")
  } else {print("failed to pass with beta distribution")}

  #test chisq distribution
  print("testing chi-square distribution")
  samp_chi <- ars(10000,dchi,a = 0.01,b=100)$sample
  ks_chi <- ks.test(samp_chi,rchisq(10000,3))
  if (expect_that(ks_chi$p.value,function(x){x>0.05}) == TRUE){
    print("passed with chisq distribution")
  } else {print("failed to pass with chisq distribution")}

  #test exp(x^2) distribution
  print("testing exp(x^2) distribution")
  if(inherits(try(ars(10000,dx2,a=-10,b=10)), "try-error")) {
    print(paste("fails to evaluate with exp(x^2) distribution"))
  }
})
