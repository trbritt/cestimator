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