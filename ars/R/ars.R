#' Adaptive Rejection Sampling
#'
#' @param sample_n A number, denoting the number of sample points generated.
#' @param g Univariate log-concave probability density function f(x).
#' @param a First starting point where f(x) is defined, default value is -Inf.
#' @param b Second starting point where f(x) is defined, default value is Inf.
#'
#' @return A list containing Information about:
#' \item{f}{logf(x)}
#' \item{sample}{A vector of size N containing sample points}
#' @author Jue Wang, Li Yang, Winnie Gao
#' @references Gilks, W.R., P. Wild. (1992) Adaptive Rejection Sampling for Gibbs Sampling, Applied Statistics 41:337–348.
#' @examples
#' #Sample from Normal Distribution
#' library(ars)
#' #Normal Distribution Density Function with mean=0 and variance=1
#' dn <- function(x){
#'    return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
#' }
#' samp_x <- ars(10000, dn, -10, 10)$sample
#' x <- seq(-10,10,len=100)
#' y <- dnorm(x)
#' plot(density(samp_x), type="l", col="red", main="Adaptive Rejection Sampling from Normal Distribution")
#' points(x,y)
#'

ars <- function(samp_n, g, a=-Inf, b=Inf){
  ##Check if the input is valid or not
  if(class(g) != "function"){
    stop('g should be a function')
  }

  require(assertthat)
  if(!assert_that(is.numeric(a))){
    stop('please provide valid lower bound')
  }

  if(!assert_that(is.numeric(b))){
    stop('please provide valid upper bound')
  }

  if(!assert_that(is.numeric(samp_n))){
    stop('please provide valid sample size')
  }

  ##Take log of input function ##
  h <- function(x){
    return(log(g(x)))
  }


  ##Initialize return value
  final_x <- rep(NA, samp_n)
  ##Set a counter to control the sample size
  counter <- 0

  ##Check starting point
  ct <- 0
  step_width <- 0.5
  start_val <- 0

  ##Case 1: a is not -Inf, b is not Inf
  if(a != -Inf && b != Inf){
    mat <- init_mat(a, b, h, Deriv_sub)
    z <- init_z(mat)
    if((abs(mat[1,3])<1e-9) && (abs(mat[2,3])<1e-9)){
      return(list(f=h, sample = runif(samp_n,a,b)))
    }
  }
  ##Case 2: a is -Inf, b is not Inf
  if(a == -Inf && b != Inf){
    if(start_val > b){
      start_val <- b - step_width
    }
    a_star <- start_val
    test_de <- Deriv_sub(a_star,fun=h,a,b)
    while( -Inf < test_de && test_de <= 0 && ct <= 50){
      a_star <- a_star - step_width
      test_de <- Deriv_sub(a_star,fun=h,a,b)
      ct <- ct + 1
    }
    mat <- init_mat(a_star, b, h, Deriv_sub)
    z <- init_z(mat) ## wrong? a_star or a?
  }
  ##Case 3: a is not -Inf, b is Inf
  if(a != -Inf && b == Inf){
    if(start_val < a){
      start_val <- a + step_width
    }
    b_star <- start_val
    test_de <- Deriv_sub(b_star, fun = h, a, b)
    while( 0 <= test_de && test_de < Inf && ct <= 50){
      b_star <- b_star + step_width
      test_de <- Deriv_sub(b_star, fun = h, a, b)
      ct <- ct + 1
    }
    mat <- init_mat(a, b_star, h, Deriv_sub)
    z <- init_z(mat)
  }
  ##Case 4: a is -Inf, b is Inf
  if(a == -Inf && b == Inf){
    a_star <- start_val - step_width
    b_star <- start_val + step_width
    test_de_a <- Deriv_sub(a_star, fun = h, a, b)
    test_de_b <- Deriv_sub(b_star, fun = h, a, b)
    while(-Inf < test_de_a && test_de_a <= 0 && ct <= 50){
      a_star <- a_star - step_width
      test_de_a <- Deriv_sub(a_star, fun = h, a, b)
      ct <- ct + 1
    }
    while(0 <= test_de_b && test_de_b < Inf && ct <= 100){
      b_star <- b_star + step_width
      test_de_b <- Deriv_sub(b_star, fun = h, a, b)
      ct <- ct + 1
    }
    mat <- init_mat(a_star, b_star, h, Deriv_sub)
    z <- init_z(mat)
  }

  ##Start Rejection Sampling process
  while(counter < samp_n){
    ##Check the log concavity of the function
    if (!check_concave(mat)){
      stop('Input must be a log-concave function')
    }
    ##Generate one sample x for each time
    samp_x_n = 1
    samp_x <- sample_x(z, mat, samp_x_n)

    ##Generate one uniformly distributed number
    w <- runif(samp_x_n)

    ##Calulate the upper bound and lowerbond at specified sample point
    u <- upper_bound(z, mat, samp_x)
    l <- lower_bound(mat, samp_x)

    ##Reject/Accept sample
    if(w <= exp(l-u)){
      counter = counter + 1
      final_x[counter] = samp_x
    }
    else{
      if (w <= exp(h(samp_x)-u)){
        counter = counter + 1
        final_x[counter] = samp_x
      }
      mat <- update_mat(mat, samp_x, h, Deriv_sub,a,b)
      z <- update_z(z, samp_x, mat)
    }
  }
  return(list(f=h, sample=final_x))
}
