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

Generally speaking, the central limit theorem does not help you much when the number of observations is short, which makes non-parametric estimation quite unideal. Benchmark estimators, like OLS, contribute to overall error via inefficiency, not usually through bias. In order to remedy this, we combine these estimators with ones that are very efficient, but display a large bias, namely constant estimators. In this way, we arrive at the James-Stein shrinkage estimator of location:

$$\widehat{\mathbf{\mu}}^S = (1-\alpha)\widehat{\mathbf{\mu}} + \alpha \mathbf{b}$$

where $\mathbf{b}$ is _any_ fixed N-dim vector (a constant estimator of location), and the optimal choice for $\alpha$ can be shown to be:

$$\alpha = \frac{1}{T}\frac{N\bar{\lambda}-2\lambda_1}{(\mathbf{\mu}-\mathbf{b})^\prime(\mathbf{\mu}-\mathbf{b})}$$

where $\lambda_1$ is the largest of $N$ eigenvalues of the $\mathbf{\Sigma}$ and $\bar{\lambda}$ is the average of the eigenvalues. The estimator can be shrunk towards the grand mean, for example, with $\mathbf{b}\to \widehat{\mathbf{\mu}}$, or the volatility-weighted grand mean, but this is part of the definition of the shrinkage you use. We use the grand-mean shrinkage target herein.

The shrinakge estimator of dispersion can be shown to be:

$$\widehat{\mathbf{\Sigma}}^S = (1-\alpha)\widehat{\mathbf{\Sigma}} + \alpha \widehat{\mathbf{C}}$$

where $\widehat{\mathbf{C}}=\frac{1}{N}\Sigma^N_{n=1}\lambda_n \mathbb{1}$ and the optimal choice of weight is:

$$\alpha = \frac{1}{T}\frac{\frac{1}{T}\Sigma^T_{t=1} \text{tr}\left[\left(\mathbf{x}_t\mathbf{x}_t^\prime - \widehat{\mathbf{\Sigma}}\right)^2\right]}{\text{tr}\left[\left(\widehat{\mathbf{\Sigma}}-\widehat{\mathbf{C}}\right)^2\right]} $$

### Maximum Likelihood

This method assumes a probability distribution of the data, and finds the maximum of the likelihood function which results in the most 
probable data. For instance, OLS estimation for linear regression maximizes likelihood when errors are normally distributed.

Let the mode of a distribution be the value that yields the peak of the distribution. Suppose that only only one observation $\mathbf{x}_1$ is available. It's reasonable to assume that it's near the mode. Therefore, after assuming a specific parametric family for the distribution, the intuitive value of the parameters are those that makes the pdf at that point the largest!

In order to apply this, we assume that the distribution is elliptically distributed with a generator $g$ and therefore the pdf will be of the form:

$$ f_\mathbf{\theta}(\mathbf{x}) =\frac{1}{\sqrt{|\mathbf{\Sigma}|}}g(\text{Ma}^2(\mathbf{x}, \mathbf{\mu}, \mathbf{\Sigma}))$$ 

where Ma is the Mahalanobis distance. You can show that the estimators obey a coupled set of implicit equations:

$$
\widehat{\mathbf{\mu}} = \sum_{t=1}^T\frac{w(\text{Ma}^2(\mathbf{x}t, \mathbf{\mu}, \mathbf{\Sigma}))}{\Sigma_{s=1}^T w(\text{Ma}^2(\mathbf{x}_s, \mathbf{\mu}, \mathbf{\Sigma}))}\mathbf{x}_t
$$

$$
  \widehat{\mathbf{\Sigma}} =\frac{1}{T}\sum_{t=1}^T (\mathbf{x}_t-\widehat{\mathbf{\mu}})^\prime w(\text{Ma}^2(\mathbf{x}_t, \mathbf{\mu}, \mathbf{\Sigma}))
$$

where the w's are defined in terms of the generator $w(z) = -2g^\prime(z)/g(z)$. Notice that these w's are weights $w_t = w(\text{Ma}^2(\mathbf{x}_t, \mathbf{\mu}, \mathbf{\Sigma}))$ showing the contribution of that observation to the resulting estimation.

### Robust Estimation

## Compilation

To build the program requires only the `Eigen3` library as a prerequisite. Other than that, just run `make` and you'll be on your way.
