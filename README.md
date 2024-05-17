# Cestimator

This is a first pass at my implementation of estimators of distributions. This will highlight a few different types of estimation.

### Non-parametric

This estimation relies solely on the raw statistics of the data, namely that the correct location of the distribution is given simply
by the mean, with the scatter matrix describing the data being the covariance of the data.

### Shrinkage

Shrinkage refers to the fact that a model describes newer data less well than the data on which it is developed, resulting in the value of
the coefficient of determination reducing. A shrinkage estimator is one that incorporates this effect, and expands on the naïve approach
of a non-parametric guess by adding some additional information to the fitting procedure.

### Maximum Likelihood

This method assumes a probability distribution of the data, and finds the maximum of the likelihood function which results in the most 
probable data. For instance, OLS estimation for linear regression maximizes likelihood when errors are normally distributed.

### Robust Estimation

## Compilation

To build the program requires only the `Eigen3` library as a prerequisite. Other than that, just run `make` and you'll be on your way.
