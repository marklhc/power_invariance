# Computing Power for LRT to Detect Measurement Non-Invariance 
2015 Mark Lai (mark.lai@uc.edu), University of Cincinnati, OH, USA

Use at your own risk

This is an R script file to compute power for likelihood ratio test to detect 
metric and scalar non-invariance with a reference and a focal group

Last modified: August 29, 2015

Required packages: lavaan, MASS

How to use:
1. Input parameters under the `User specification` section
2. Run the whole script in R

Assumptions:
+ All observed variables are standardized (variance = 1)
+ The items are unidimensional (one factor model)
+ The latent factor has unit variance for both groups
+ Multivariate normality
