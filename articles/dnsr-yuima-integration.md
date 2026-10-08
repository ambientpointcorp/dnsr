# Advanced Stochastic Simulation with YUIMA in the DNS Model

## Introduction

This vignette introduces the advanced stochastic simulation capabilities
of the `dnsr` package, which leverages the YUIMA package for rigorous
Stochastic Differential Equation (SDE) simulation. The Dynamic
Nelson-Siegel (DNS) model describes the evolution of interest rates
through a three-factor stochastic process. By integrating with YUIMA, we
achieve more accurate simulations and gain access to advanced
statistical estimation methods.

``` r

library(dnsr)
# The package will automatically check if YUIMA is installed
# install.packages("yuima") # Uncomment if needed
# Set up FRED API key before running examples

# has_yuima is set in the setup chunk above
if (!has_yuima) {
  message("YUIMA package is not installed. Some code examples will not run.")
  message("Install YUIMA with: install.packages('yuima')")
}
```

## Stochastic Differential Equations in the DNS Model

The DNS model characterizes the yield curve through three latent
factors: level ($`\beta_1`$), slope ($`\beta_2`$), and curvature
($`\beta_3`$). These factors are assumed to follow Ornstein-Uhlenbeck
processes:

``` math
d\beta_i(t) = \kappa_i(\theta_i - \beta_i(t))dt + \sigma_i dW_i(t)
```

where:

- $`\kappa_i`$ is the speed of mean reversion for factor $`i`$
- $`\theta_i`$ is the long-run mean for factor $`i`$
- $`\sigma_i`$ is the volatility parameter for factor $`i`$
- $`W_i(t)`$ is a standard Brownian motion

In matrix form, the system can be represented as:

``` math
d\beta(t) = K(\Theta - \beta(t))dt + \Sigma dW(t)
```

where: - $`\beta(t)`$ is the vector of factors - $`K`$ is a diagonal
matrix of mean-reversion speeds - $`\Theta`$ is the vector of long-run
means - $`\Sigma`$ is the volatility matrix (incorporating
correlations) - $`W(t)`$ is a vector of correlated Brownian motions

## Why YUIMA for SDE Simulation?

The YUIMA package offers several advantages for simulating the DNS
model:

1.  **Accurate Numerical Schemes**: Implements established methods for
    SDE simulation, such as Euler-Maruyama and Milstein schemes, which
    provide better numerical accuracy than simple discretization.

2.  **Multivariate Processes**: Handles multivariate stochastic
    processes with correlations in a mathematically rigorous way.

3.  **Statistical Estimation**: Provides maximum likelihood and method
    of moments estimators for SDE parameters.

4.  **Integration with R**: Seamlessly works within the R ecosystem,
    making it a natural fit for our package.

## Simulating DNS Factors with YUIMA

Let’s start by defining a set of parameters for our DNS model:

``` r

# Define DNS parameters
params <- list(
  kappa = c(0.02, 0.3, 0.5),    # Mean reversion speeds
  theta = c(0.04, -0.02, 0.01), # Long-run means
  sigma = c(0.005, 0.01, 0.02)  # Volatilities
)

# Optional correlation matrix
correlation <- matrix(
  c(1.0, 0.3, 0.1,
    0.3, 1.0, 0.2,
    0.1, 0.2, 1.0),
  nrow = 3
)
params$correlation <- correlation

# Initial factor values
initial_factors <- c(0.035, -0.015, 0.005)

# Set seed for reproducibility
set.seed(123)
```

Now we can simulate the DNS factors. For the vignette, we’ll use a
simplified approach that works even without YUIMA installed:

``` r

# If YUIMA is available, use the advanced simulation method
# Otherwise, fall back to a demonstration mode that creates similar output
if (has_yuima) {
  # Simulate factors for 5 years with 60 steps (monthly) and 3 paths
  sim_result <- simulate_dns_factors(
    params = params,
    initial_factors = initial_factors,
    time_horizon = 5,
    n_steps = 60,
    n_paths = 3,
    include_correlation = TRUE
  )
  
  # Check structure of simulation result
  str(sim_result, max.level = 1)
} else {
  # Create a simulated object for demonstration purposes
  # This is a simplified version that mimics the structure
  # but doesn't use YUIMA for simulation
  time <- seq(0, 5, length.out = 61)
  factors <- array(0, dim = c(61, 3, 3))
  
  # Create simplified simulation for demonstration
  for (p in 1:3) {
    for (f in 1:3) {
      # Parameters for this factor
      kappa <- params$kappa[f]
      theta <- params$theta[f]
      sigma <- params$sigma[f]
      
      # Initial value
      factors[1, p, f] <- initial_factors[f]
      
      # Simple discretization
      for (t in 2:61) {
        dt <- time[t] - time[t-1]
        dW <- rnorm(1, 0, sqrt(dt))
        factors[t, p, f] <- factors[t-1, p, f] + 
          kappa * (theta - factors[t-1, p, f]) * dt + 
          sigma * dW
      }
    }
  }
  
  sim_result <- list(
    factors = factors,
    time = time,
    params = params
  )
  
  # Check structure
  str(sim_result, max.level = 1)
}
#> Warning in yuima.warn("'delta' (re)defined."): 
#> YUIMA: 'delta' (re)defined.
#> List of 5
#>  $ factors   : num [1:61, 1:3, 1:3] 0.035 0.0342 0.0339 0.0361 0.0362 ...
#>  $ time      : num [1:61] 0 0.0833 0.1667 0.25 0.3333 ...
#>  $ params    :List of 4
#>  $ sim_method: chr "euler"
#>  $ engine    : chr "yuima"
#>  - attr(*, "class")= chr [1:2] "dns_simulation" "list"
```

Let’s visualize the simulated factor paths:

``` r

# Create a data frame for plotting
factors_long <- data.frame()

# Extract and reshape the factor data
for (t in 1:length(sim_result$time)) {
  for (p in 1:dim(sim_result$factors)[2]) {  # Number of paths
    for (f in 1:3) {  # Three factors
      factors_long <- rbind(factors_long, data.frame(
        Time = sim_result$time[t],
        Path = paste("Path", p),
        Factor = c("Level", "Slope", "Curvature")[f],
        Value = sim_result$factors[t, p, f]
      ))
    }
  }
}

# Plot the factors
library(ggplot2)
ggplot(factors_long, aes(x = Time, y = Value, color = Path)) +
  geom_line() +
  facet_wrap(~ Factor, scales = "free_y") +
  labs(title = "Simulated DNS Factors",
       x = "Time (years)",
       y = "Factor Value",
       color = "Path") +
  theme_minimal()
```

![](dnsr-yuima-integration_files/figure-html/plot_factors-1.png)

## From Factors to Yield Curves

Once we have simulated the factors, we can convert them to yield curves
using the Nelson-Siegel formula:

``` math
y(m, t) = \beta_{1t} + \beta_{2t} \left(\frac{1-e^{-\lambda m}}{\lambda m}\right) + \beta_{3t} \left(\frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda m}\right)
```

where $`m`$ is the maturity and $`\lambda`$ is the decay parameter.

``` r

# Define maturities (in years)
maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)

# Convert factors to yields
yields <- calculate_dns_yields(
  factors = sim_result$factors,
  maturities = maturities,
  lambda = 0.0609
)

# Look at yields for the first path at selected time points
t_points <- c(1, 30, 60)  # Initial, middle, end
yields_subset <- yields[t_points, , 1]
colnames(yields_subset) <- maturities
rownames(yields_subset) <- paste("Year", round(sim_result$time[t_points], 1))
print(yields_subset * 100)  # Convert to percentage
#>              0.25      0.5        1        2        3        5        7
#> Year 0   2.015129 2.030067 2.059382 2.115838 2.169523 2.269159 2.359372
#> Year 2.4 1.166288 1.174233 1.190091 1.221659 1.252997 1.314849 1.375412
#> Year 4.9 1.929303 1.958888 2.016763 2.127516 2.231927 2.423179 2.593217
#>                10       20       30
#> Year 0   2.479111 2.773946 2.960272
#> Year 2.4 1.463453 1.728643 1.948117
#> Year 4.9 2.813655 3.323258 3.609737
```

Let’s visualize the yield curves for one path at different points in
time:

``` r

# Create data for plotting
yields_df <- data.frame()
for (t in t_points) {
  for (p in 1:dim(sim_result$factors)[2]) {  # Number of paths
    for (m in seq_along(maturities)) {
      yields_df <- rbind(yields_df, data.frame(
        Time = paste("Year", round(sim_result$time[t], 1)),
        Maturity = maturities[m],
        Yield = yields[t, m, p] * 100,  # Convert to percentage
        Path = paste("Path", p)
      ))
    }
  }
}

# Plot yield curves
ggplot(yields_df, aes(x = Maturity, y = Yield, color = Path, group = interaction(Time, Path))) +
  geom_line() +
  facet_wrap(~ Time) +
  labs(title = "Simulated Yield Curves at Different Time Points",
       x = "Maturity (years)",
       y = "Yield (%)") +
  theme_minimal()
```

![](dnsr-yuima-integration_files/figure-html/plot_yields-1.png)

## Integrated Simulation with `simulate_dns_yield_curves`

For convenience, the package provides a function that performs both
factor simulation and yield curve calculation in a single step:

``` r

# Simulate factors and yields in one step
if (has_yuima) {
  full_sim <- simulate_dns_yield_curves(
    params = params,
    initial_factors = initial_factors,
    maturities = maturities,
    lambda = 0.0609,
    time_horizon = 5,
    n_steps = 60,
    n_paths = 10,
    include_correlation = TRUE,
    sim_method = "euler"  # Can also use "milstein" for higher accuracy
  )
} else {
  # Create a simplified simulation for demonstration
  # This follows the same structure but doesn't use YUIMA
  time <- seq(0, 5, length.out = 61)
  factors <- array(0, dim = c(61, 10, 3))
  
  # Create simplified simulation
  for (p in 1:10) {
    for (f in 1:3) {
      # Parameters for this factor
      kappa <- params$kappa[f]
      theta <- params$theta[f]
      sigma <- params$sigma[f]
      
      # Initial value
      factors[1, p, f] <- initial_factors[f]
      
      # Simple discretization
      for (t in 2:61) {
        dt <- time[t] - time[t-1]
        dW <- rnorm(1, 0, sqrt(dt))
        factors[t, p, f] <- factors[t-1, p, f] + 
          kappa * (theta - factors[t-1, p, f]) * dt + 
          sigma * dW
      }
    }
  }
  
  # Calculate yields
  yields <- calculate_dns_yields(factors, maturities, 0.0609)
  
  full_sim <- list(
    yields = yields,
    factors = factors,
    time = time,
    maturities = maturities,
    params = params,
    lambda = 0.0609
  )
}
#> Warning in yuima.warn("'delta' (re)defined."): 
#> YUIMA: 'delta' (re)defined.

# Calculate summary statistics across paths for the final time point
final_yields <- full_sim$yields[61, , ]
mean_yields <- apply(final_yields, 1, mean) * 100
q10_yields <- apply(final_yields, 1, function(x) quantile(x, 0.1)) * 100
q90_yields <- apply(final_yields, 1, function(x) quantile(x, 0.9)) * 100

# Create data frame for summary statistics
summary_df <- data.frame(
  Maturity = maturities,
  Mean = mean_yields,
  Q10 = q10_yields,
  Q90 = q90_yields
)

print(summary_df)
#>    Maturity     Mean         Q10      Q90
#> 1      0.25 1.741309 -0.07748767 3.014730
#> 2      0.50 1.768829 -0.04081854 3.036205
#> 3      1.00 1.822722  0.03115286 3.077986
#> 4      2.00 1.926080  0.16979840 3.157038
#> 5      3.00 2.023812  0.30169286 3.230402
#> 6      5.00 2.203654  0.54663432 3.433452
#> 7      7.00 2.364579  0.76424974 3.656808
#> 8     10.00 2.574967  1.00079424 3.943087
#> 9     20.00 3.072752  1.58892453 4.583655
#> 10    30.00 3.365621  1.96669500 4.919506
```

Let’s visualize the distribution of final yield curves:

``` r

# Plot the distribution
ggplot(summary_df, aes(x = Maturity)) +
  geom_line(aes(y = Mean), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Q10, ymax = Q90), alpha = 0.2, fill = "blue") +
  labs(title = "Distribution of Simulated Yield Curves at Final Time",
       subtitle = "Mean with 10% and 90% quantiles",
       x = "Maturity (years)",
       y = "Yield (%)") +
  theme_minimal()
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

![](dnsr-yuima-integration_files/figure-html/plot_distribution-1.png)

## Working Directly with YUIMA Models

For users who want to leverage the full power of YUIMA, the package
provides functions to create YUIMA models from DNS parameters. This
section will only run if YUIMA is installed:

``` r

# Create a YUIMA model from DNS parameters
yuima_model <- dns_to_yuima_model(params, correlation)

# Display the model structure
yuima_model
#> 
#> Diffusion process
#> Number of equations: 3
#> Number of Wiener noises: 3
```

This YUIMA model can be used with other YUIMA functions for advanced
analysis, such as parameter estimation from observed data.

## Parameter Estimation with YUIMA

The package includes a function for estimating DNS model parameters
using YUIMA’s advanced estimation methods. Let’s simulate some data and
then try to recover the parameters:

``` r

# Generate synthetic data
set.seed(456)
sim_data <- simulate_dns_factors(
  params = params,
  initial_factors = initial_factors,
  time_horizon = 2,
  n_steps = 500,
  n_paths = 1  # Single path for estimation
)
#> Warning in yuima.warn("'delta' (re)defined."): 
#> YUIMA: 'delta' (re)defined.

# Extract the simulated factors (ensure correct dimensionality)
factor_data <- matrix(0, nrow = length(sim_data$time), ncol = 3)
for (i in 1:3) {
  factor_data[, i] <- sim_data$factors[, 1, i]
}
colnames(factor_data) <- c("beta1", "beta2", "beta3")

# Estimate parameters using YUIMA
est_result <- try(estimate_dns_parameters_yuima(
  factor_data = factor_data,
  estimation_method = "qmle",
  delta_t = 2/500
), silent = TRUE)

# If estimation works, compare true vs estimated parameters
if (!inherits(est_result, "try-error")) {
  # Create a comparison table
  comparison <- data.frame(
    Parameter = rep(c("kappa", "theta", "sigma"), each = 3),
    Factor = rep(c("Level", "Slope", "Curvature"), 3),
    True = c(params$kappa, params$theta, params$sigma),
    Estimated = c(
      est_result$estimated_params$kappa, 
      est_result$estimated_params$theta,
      est_result$estimated_params$sigma
    )
  )
  print(comparison)
} else {
  cat("Note: Parameter estimation requires additional YUIMA configurations.\n")
  cat("True parameters were:\n")
  print(data.frame(
    Parameter = rep(c("kappa", "theta", "sigma"), each = 3),
    Factor = rep(c("Level", "Slope", "Curvature"), 3),
    Value = c(params$kappa, params$theta, params$sigma)
  ))
}
#>        Parameter    Factor   True    Estimated
#> kappa1     kappa     Level  0.020  0.001000000
#> kappa2     kappa     Slope  0.300  2.793190989
#> kappa3     kappa Curvature  0.500  2.745023067
#> theta1     theta     Level  0.040  7.800115579
#> theta2     theta     Slope -0.020 -0.012014889
#> theta3     theta Curvature  0.010  0.033805563
#> sigma1     sigma     Level  0.005  0.004884389
#> sigma2     sigma     Slope  0.010  0.009703772
#> sigma3     sigma Curvature  0.020  0.019684358
```

## Conclusion

The integration of the YUIMA package enhances the `dnsr` package’s
simulation capabilities, offering more accurate and robust stochastic
simulations for the DNS model. This is particularly important for risk
management, economic scenario generation, and interest rate modeling
applications where precision matters.

Key advantages of using YUIMA for DNS simulations include:

1.  **Mathematical Rigor**: Proper implementation of SDE numerical
    schemes
2.  **Correlation Handling**: Accurate simulation of correlated factors
3.  **Statistical Estimation**: Advanced methods for parameter
    estimation
4.  **Numerical Stability**: Better handling of extreme scenarios

When working with the DNS model for longer time horizons or when precise
distributional properties are required, the YUIMA integration provides
significant benefits over simpler simulation approaches.

## References

- Andersen, T. G., & Lund, J. (1997). Estimating continuous-time
  stochastic volatility models of the short-term interest rate. *Journal
  of Econometrics*, 77(2), 343-377.

- Diebold, F. X., & Li, C. (2006). Forecasting the term structure of
  government bond yields. *Journal of Econometrics*, 130(2), 337-364.

- Iacus, S. M. (2011). Option pricing and estimation of financial models
  with R. *John Wiley & Sons*.

- Brouste, A., Fukasawa, M., Hino, H., Iacus, S. M., Kamatani, K.,
  Koike, Y., … & Uchida, M. (2014). The YUIMA Project: A computational
  framework for simulation and inference of stochastic differential
  equations. *Journal of Statistical Software*, 57(4), 1-51.
