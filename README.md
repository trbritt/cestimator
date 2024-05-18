# Cestimator

This is a first pass at my implementation of estimators of distributions. This will highlight a few different types of estimation.

### Non-parametric

This estimation relies solely on the raw statistics of the data. What I mean by that is that we're essentially trying to follow the _law of large numbers_, where:

$$\lim_{T\to\infty}\frac{1}{T}\sum_{t=1}^T [ \text{past} ] \approx \mathbb{E}[\text{future}] $$

This directly implies the Glivenko-Cantelli theorem, which states that the distribution of a set of iid variables tends to the true distribution in the large observation limit. So for example, in terms of the cumulative distribution functions:

$$\lim_{T\to\infty} F_{i_T}(\mathbf{x}) = F_\mathbf{X}(\mathbf{x})$$

where $i_T$ is the set of previous realisations of the variable of interest for a history period T. This brings us to say that the distribution of the previous information is simply the sum of the previous realisiations of the variable, since those outcomes are definitely possible:

$$f_{i_T}(\mathbf{x}) = \frac{1}{T} \sum_{t=1}^T \delta^{\mathbf{x}_t}(\mathbf{x})$$

This lets us derive the non-parametric estimator of the expected value, producing the sample mean:

$$\widehat{E}[i_T] = \sum_{\mathbb{R}^N}\mathbf{x}f_{i_T}(\mathbf{x})d\mathbf{x} = \frac{1}{T}\sum_{t=1}^T\mathbf{x}_t$$

Likewise, the non-parametric estimator of dispersion:

$$\widehat{Cov}[i_T] = \frac{1}{T}\sum_{t=1}^T \left(\mathbf{x}_t - \widehat{E}[i_T])\right)\left(\mathbf{x}_t - \widehat{E}[i_T])\right)^\prime$$

with transposition denoted by $\prime$. The geometric interpretation of these parameters is that they define an ellipse through the multidimensional space where the average Mahalnobis, for the given location and dispersion, be distance of 1.

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
