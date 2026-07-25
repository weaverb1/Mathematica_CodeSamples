# In this problem, you will analyze data from a Weibull distribution by 
# estimating parameters and evaluating likelihoods using numerical methods.

library(tidyverse)

# Part a: Weibull Distribution Function

# Write a function that uses numerical integration to calculate the expected 
# value and variance of a Weibull distribution. The function should accept the 
# shape (k) and scale (λ) parameters as inputs.

# This weibull function handles the issue with x's upper bound being infinity by 
# choosing a sufficiently large b to act as the upper bound. When b and n are 
# chosen well, it works great, but I wanted to find another method that was less 
# dependent on choosing a good bound and sample size which you will see in the 
# second function I created right after. 
weibull <- function(k, lambda, n = 100, b = 1000){
  a <- 0
  h <- (b-a)/n
  
  pdf <- function(x, k, lambda) (k/lambda)*((x/lambda)^(k-1))*exp(-(x/lambda)^k)
  e_x <- function(x, k, lambda) x*pdf(x,k,lambda)
  e_x2 <- function(x, k, lambda) (x^2)*pdf(x,k,lambda)
  
  x_vals <- seq(a, b, length = n + 1)
  e_x_vals <- e_x(x_vals, k, lambda)
  e_x2_vals <- e_x2(x_vals, k, lambda)
  
  mean <- h / 2 * (e_x_vals[1] + 2 * sum(e_x_vals[2:n]) + e_x_vals[n + 1])
  e_x2_final <- h / 2 * (e_x2_vals[1] + 2 * sum(e_x2_vals[2:n]) + e_x2_vals[n + 1])
  var <- e_x2_final - (mean)^2
  
  return(c(mean,var))
}

#code to check the accuracy of my function
k <- 1
lambda <- 3

true_mean <- lambda*gamma(1 + 1/k)
true_variance <- (lambda^2) * (gamma(1 + 2 / k) - (gamma(1 + 1 / k))^2)

weibull(k, lambda, b =100)
true_mean
true_variance

# With help from chatGPT I was able to learn how to map the weibull function to 
# a finite interval of 0 to 1. This function works well and does not require the 
# user to pick a smart upper bound to get a reliable estimate.
weibull <- function(k, lambda, n = 100){
  
  
  pdf <- function(x, k, lambda) (k/lambda)*((x/lambda)^(k-1))*exp(-(x/lambda)^k)
  
  g_ex <- function(t, k, lambda){
    x <- t/(1-t)
    x*(pdf(x,k,lambda)/(1-t)^2)
  }
  
  g_ex2 <- function(t, k, lambda){
    x <- t/(1-t)
    (x^2)*(pdf(x,k,lambda)/(1-t)^2)
  }
  
  h <- 1/n
  t <- seq(0, 1, length.out = n + 1)[-n - 1] 
  
  e_x_vals <- sapply(t, g_ex, k = k, lambda = lambda)
  e_x2_vals <- sapply(t, g_ex2, k = k, lambda = lambda)
  
  mean <- h * (sum(e_x_vals) - 0.5 * (e_x_vals[1] + e_x_vals[length(e_x_vals)]))
  e_x2_final <- h * (sum(e_x2_vals) - 0.5 * (e_x2_vals[1] + e_x2_vals[length(e_x2_vals)]))
  var <- e_x2_final - (mean)^2
  
  return(c(mean,var))
}

#code to check the accuracy of my function
k <- 1
lambda <- 3

true_mean <- lambda*gamma(1 + 1/k)
true_variance <- (lambda^2) * (gamma(1 + 2 / k) - (gamma(1 + 1 / k))^2)

weibull(k, lambda)
true_mean
true_variance

# Part b: Parameter Estimation using Newton-Raphson

# Use the data in problem2dat.csv where the only variable is assumed to follow a
# Weibull distribution with a fixed shape parameter k.

# Write a function that uses the Newton-Raphson method to maximize the 
# log-likelihood of the Weibull distribution by optimizing only the scale 
# parameter λ.

# The function should return both the log-likelihood and the estimated scale 
# parameter λ.

# this function creates the log-likelihood and its derivatives needed for the MLE 
# function
make_functions <- function(data, k){
  n <- length(data)
  r <- sum(log(data))
  lr <- sum(data^k)
  
  log.likelihood <- function(lambda){
    n*log(k) + (k-1)*r -n*k*log(lambda) - (lambda^(-k))*lr
  }
  
  d1_log_like <- function(lambda){
    -((n*k)/lambda) + k*((lambda)^(-k-1))*lr
  }
  
  d2_log_like <- function(lambda){
    (n*k)/(lambda^2)-k*(k+1)*(lambda^(-k-2))*lr
  }
  list(logLike = log.likelihood, d1 = d1_log_like, d2 = d2_log_like)
}

# This function uses Newton-Rhapson to optimize the variable lambda for a given 
# k for a weibull distribution
mle_weibull <- function(data, lambda, k, tolerance = 1e-06, max_iter = 1000){
  if (lambda < 0) stop("Error: Starting lambda must be greater than 0")
  if (k < 0) stop("Error: Value 'k' must be greater than 0")
  if (min(data) < 0) stop("Error: Data values must be 0 or greater")
  
  f <- make_functions(data, k)
  count <- 0
  
  while (abs(f$d1(lambda)) > tolerance && count < max_iter){
    if (is.nan(f$d1(lambda)) || is.nan(f$d2(lambda))) stop("Error: Numerical instability encountered.")
    lambda <- lambda - f$d1(lambda)/f$d2(lambda)
    if (lambda <= 0) stop("Error: lambda became non-positive during iterations.")
    count <- count + 1
  }
  if (count >= max_iter) warning("MLE did not converge within the iteration limit")
  if (f$d2(lambda) >= 0) warning("MLE did not maximize. Try a different starting lambda.")
  list(estimated_lambda = lambda, log_likelihood = f$logLike(lambda))
}

# ckeck function
check <- rweibull(100, shape = 2, scale = 5)
k <- 2

# code suggested by chat while debugging to calculate a thoughtful starting value
lambda_init <- mean(check) / gamma(1 + 1 / k)
result <- mle_weibull(check, lambda_init, k)
print(result)

# using provided data
read_data <- read.csv("problem2dat.csv")
data <- read_data$V
k <- 2
lambda_init <- mean(data) / gamma(1 + 1 / k)

result <- mle_weibull(data, lambda_init, k)
print(result)

# Part c: Grid Search for Parameter Estimation

# Conduct a grid search over the shape parameter k, varying it from 1 to 4 in 
# increments of 0.1.

# For each value of k in the grid, use the Newton-Raphson method to find the 
# optimal scale parameter λ that maximizes the log-likelihood for that fixed k.

# Identify and report the (k,λ) pair that maximizes the log-likelihood.

# this code finds the optimized lambda for each k along with the log-likelihood
k <- seq(1,4, by = 0.1)
lambda_init <- mean(data) / gamma(1 + 1 / k)
results <- mapply(function(k, lambda_init) mle_weibull(data, lambda_init,k), k, lambda_init, SIMPLIFY = FALSE)

# this code creates a data frame with the different parts of our grid search
results_df <- data.frame(
  k = k,
  lambda = sapply(results, function(res) res$estimated_lambda),
  log_likelihood = sapply(results, function(res) res$log_likelihood)
)

# this code finds the (k, lambda) pair with the highest log-likelihood
max_k_lam <- results_df |>
  filter(log_likelihood == max(log_likelihood))

print(max_k_lam)

paste0("The (k,λ) pair that maximizes the log-likelihood is: ", round(max_k_lam$k,3), ", ", round(max_k_lam$lambda,3))
