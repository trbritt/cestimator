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