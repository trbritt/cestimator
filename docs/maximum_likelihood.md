
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
