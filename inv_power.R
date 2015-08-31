# ----------------------------------------------------------------------- #
# 2015 Mark Lai (mark.lai@uc.edu), University of Cincinnati, OH, USA
#
# Use at your own risk
# 
# Script to compute power for likelihood ratio test to detect 
# metric and scalar non-invariance with a reference and a focal group
#
# Last modified: August 29, 2015
#
# Required packages: lavaan, MASS
# 
# How to use:
#   1. Input parameters under the `User specification` section
#   2. Run the whole script in R
#
# Assumptions:
#   - All observed variables are standardized (variance = 1)
#   - The items are unidimensional (one factor model)
#   - The latent factor has unit variance for both groups
#   - Multivariate normality
# ----------------------------------------------------------------------- #


# User specification ------------------------------------------------------

p        <- 6    # number of items
pninv    <- 1    # number of non-invariant items
dif      <- 0.5  # magnitude of non-invariance; can be interpreted as Cohen's d
lambda   <- .7   # (average) factor loading / pattern coefficient
d_kappa  <- 0.5  # difference in factor means; can be interpreted as Cohen's d
N_r      <- 100  # sample size for the reference group
N_f      <- 100  # sample size for the focal group
alpha    <- .05  # Significance level for LRT
MCsample <- 500  # number of Monte Carlo samples (recommended 500 or more)

# End User specification --------------------------------------------------

# Generate population covariances for the reference and the focal group
phi <- 1  # factor variance (assumed constant across groups)
Lambda_r <- rep(lambda, p)  # Lambda matrix for the reference group
Lambda_f <- c(rep(lambda, p - pninv), 
              rep(lambda - dif, pninv))  # Lambda matrix for the focal group
Sigma_r <- Lambda_r %*% as.matrix(phi) %*% t(Lambda_r)  # population covariance
diag(Sigma_r) <- 1                                      # of reference group
Sigma_f <- Lambda_f %*% as.matrix(phi) %*% t(Lambda_f)  # population covariance
diag(Sigma_f) <- 1                                      # of focal group

# Generate sample covariances
## Sample covariance matrix for reference group
varnames <- paste0("x", seq_len(p))
S_r <- rWishart(MCsample, df = N_r, Sigma = Sigma_r) / (N_r - 1)
## Sample covariance matrix for focal group
S_f <- rWishart(MCsample, df = N_f, Sigma = Sigma_f) / (N_f - 1)
colnames(S_r) <- colnames(S_f) <- varnames

# Define metric invariance model
cfa_model <- 
  paste('f1 =~', paste(varnames, collapse = " + "))
library(lavaan)
power_vec_m <- numeric(MCsample)
for (i in seq_len(MCsample)) {
  config_fit <- cfa(cfa_model, sample.cov = list(S_r[ , , i], S_f[ , , i]), 
                    sample.nobs = list(N_r, N_f))
  metric_fit <- update(config_fit, group.equal = c("loadings"))
  power_vec_m[i] <- anova(config_fit, metric_fit)$"Pr(>Chisq)"[2]
}

metric_power <- mean(power_vec_m < alpha)

# Generate mean vectors for the reference and the focal group
kappa_r <- 0  # factor mean (assumed constant across groups)
kappa_f <- kappa_r - d_kappa
tau <- 0  # measurement intercept
tau_r <- rep(tau, p)  # tau vector for the reference group
tau_f <- c(rep(tau, p - pninv), 
           rep(tau - dif, pninv))  # tau vector for the reference group
mu_r <- tau_r + Lambda_r %*% as.matrix(kappa_r)
mu_f <- tau_f + Lambda_r %*% as.matrix(kappa_f)

# Generate sample covariance for focal group (with metric invariance)
S_f <- rWishart(MCsample, df = N_f, Sigma = Sigma_r) / (N_r - 1)
colnames(S_f) <- varnames

# Generate sample mean vectors for both groups
M_r <- MASS::mvrnorm(MCsample, mu_r, Sigma = Sigma_r / N_r)
M_f <- MASS::mvrnorm(MCsample, mu_f, Sigma = Sigma_r / N_f)

power_vec_s <- numeric(MCsample)
for (i in seq_len(MCsample)) {
  metric_fit <- cfa(cfa_model, meanstructure = TRUE, 
                    sample.cov = list(S_r[ , , i], S_f[ , , i]), 
                    sample.mean = list(M_r[i, ], M_f[i, ]), 
                    sample.nobs = list(N_r, N_f), group.equal = c("loadings"))
  scalar_fit <- update(metric_fit, sample.cov = list(S_r[ , , i], S_f[ , , i]), 
                       sample.mean = list(M_r[i, ], M_f[i, ]), 
                       group.equal = c("loadings", "intercepts"))
  power_vec_s[i] <- anova(metric_fit, scalar_fit)$"Pr(>Chisq)"[2]
}

scalar_power <- mean(power_vec_s < alpha)

cat("\n=======================================================================", 
    "\n| Result:\n| ", 
    "\n| Power for detecting", pninv, "out of", p, 
    "loading difference(s) of magnitude", dif, "=", metric_power, "\n| ", 
    "\n| Power for detecting", pninv, "out of", p, 
    "intercept difference(s) of magnitude", dif, "=", scalar_power, 
    "\n=======================================================================\n")