# In this problem, you’ll examine the impact of collinearity on the F statistic 
# in a simple linear regression model with two predictors.

# libraries
library(MASS)
library(tidyverse)

# Part a: Generating Bivariate Normal Data

# Write a function to generate a dataset with n = 40 samples and two predictors,
# X1 and X2, where the correlation between X1 and X2 is specified by a parameter
# ρ. The response variable Y should be generated based on the relationship 
# Y = 0.5X1 − 0.5X2 + ϵ, where ϵ ∼ N(0,1) is a normally distributed error term.

generate_data <- function(rho, n = 40) {
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  X <- mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  error <- rnorm(n)
  Y <- 0.5*X[,1] - 0.5*X[,2] + error
  return(data.frame(Y = Y, X1 = X[,1], X2 = X[,2]))
}

#checking to make sure the input is in the expected form
check <- generate_data(.8, 40)
check

# Part b: F Statistic Calculation Function

# Write a function to compute the F statistic for a linear regression model, 
# given a dataset with predictors X1, X2, and response Y .

# using lm to calculate y_hat
F_stat <- function(data){
  Y <- data$Y
  X1 <- data$X1
  X2 <- data$X2
  mod <- lm(Y ~ X1 + X2, data = data)
  y_hat <- predict(mod)
  y_bar <- mean(Y)
  
  ssr <- sum( (y_hat - y_bar)^2 )
  sse <- sum( (Y - y_hat)^2 )
  
  f_stat <- (ssr/2)/(sse/(length(Y)-2-1))
  return(f_stat)
}

# Compare results between my code and pre-built r functions
f_check <- F_stat(check)
mod <- lm(Y ~ X1 + X2, data = check)
f_summary <- summary(mod)$fstatistic[1]

# Display
f_check
f_summary

# Part c: Bootstrapping

# Simulate a data set with ρ = 0.8, use your function to calculate the F statistic.
# Derive a distribution of the F statistic using bootstrapping on that data set
# Report the 95% confidence interval of the F statistic based on your bootstrapping results.

bootstrap_ci <- function(rho = 0.8, iter = 1000, n = 40){
  data <- generate_data(rho, n=n)
  f_stats <- numeric(iter)
  for(i in 1:iter){
    samp_data <- data[sample(nrow(data), replace = TRUE), ]
    f_stats[i] <- F_stat(samp_data)
  }
  ci <- quantile(f_stats, c(0.025, 0.975))
  return(ci)
}

ci <- bootstrap_ci(0.8, 100)
paste0("Bootstrap Conf. Int: ", round(ci[1],3),", ",round(ci[2],3))

# Part d: Simulation Study

# Conduct a simulation study where you vary the correlation ρ from 0 to 0.9 in 
# increments of 0.1, generating a new dataset for each value of ρ and 
# calculating the F statistic each time.

# Analyze how the F statistic changes with increasing collinearity.

# this function calculates n f-statistics for different datasets generated at 
# each rho value and returns the f-stat with its rho value.
sim <- function(n = 5){
  rho_values <- seq(0, 0.9, by = 0.1)
  rho <- rep(rho_values, each = n)
  f <- numeric(length(rho))
  for(i in seq_along(rho)){
    data <- generate_data(rho[i])
    f[i] <- F_stat(data)
  }
  final_data <- data.frame(rho = rho, f_stat = f) |>
    arrange(rho)
  return(final_data)
}

# I created plots with different n values to see what general patterns the 
# f-stat followed for different rho values.
sim_dat <- sim(10)
plot(sim_dat$rho, sim_dat$f_stat, main = "How rho effects the f-statistic", xlab= "rho", ylab="F-statistic")

sim_dat <- sim(50)
plot(sim_dat$rho, sim_dat$f_stat, main = "How rho effects the f-statistic", xlab= "rho", ylab="F-statistic")

sim_dat <- sim(100)
plot(sim_dat$rho, sim_dat$f_stat, main = "How rho effects the f-statistic", xlab= "rho", ylab="F-statistic")

# The results of the simulation study suggest that as rho increases the range of
# possible f-statistics decreases. We can see this in the above plots. In the 
# plots we can see that the y-values for rho = 0 have a much higher maximum that 
# rho = 0.9. This suggests that models that include multiple highly correlated 
# variables will tend to do worse at explaining the variation in the dependent 
# variable than models with less correlated variables. This fits with the 
# instructions that we have received from Dr. Fisher in Stat 536, where he has 
# warned us to check for and avoid including correlated variables in our models 
# as it will often weaken the model.

# Part e: External Data Analysis

ex_data <- read.csv("problem1dat.csv")

# Compute the F statistic for this dataset.

paste0("F-stat for the provided dataset: ", round(F_stat(ex_data),3))

# Part f: Randomization Test:

# Perform a randomization test to test the hypothesis that all the regression 
# coefficients β1 and β2 are 0, using the F statistic as your test statistic.

# In the randomization test, shuffle the response variable Y multiple times 
# (e.g., 1000 times), recomputing the F statistic for each shuffled dataset.

# Use the distribution of the F statistic under randomization to assess whether 
# the observed F statistic for the original data is statistically significant.

set.seed(1312)
random_test <- function(data, iter = 1000){
  observed_f <- F_stat(data)
  ran_f <- numeric(length(iter))
  
  for(i in 1:iter){
    y_shuffle <- sample(data$Y, replace = FALSE)
    shuffled_data <- data.frame(Y = y_shuffle, X1 = data$X1, X2 = data$X2)
    ran_f[i] <- F_stat(shuffled_data)
  }
  p_val <- mean(ran_f >= observed_f)
  return(p_val)
}

paste0("With a p-value of ", random_test(ex_data), " we can determine that the f-stat for the original data is statistically significant.")
