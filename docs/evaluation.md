## Evaluation

So far, the code has consisted of multiple estimators through which one can analayse the relationship between potentially multiple sources of bias and inefficiency. The goal, now is to take these analyses for many potential sources, and to combine them to optimize a given objective. For example, the value of a portfolio $w_T = \boldsymbol{\alpha}^\prime\mathbf{p}_T$ is the made at decision time $T$ of prices $p_T^{(n)},n\in[0,n]$, each price being allocated a weight $\alpha_n$. The projection of the portfolio to an investment horizon $\tau$ is therefore a random variable $W_{T+\tau}(\boldsymbol{\alpha}) = \boldsymbol{\alpha}^\prime \mathbf{P}_{T+\tau}$.

In all cases, the objective functional $\Psi_{\boldsymbol{\alpha}}=\boldsymbol{\alpha}^\prime\mathbf{M}$ will be a linear function of the market vector $\mathbf{M}=\mathbf{a}+\mathbf{B}\mathbf{P}_{T+\tau}$. 

- For $\mathbf{a}=\mathbf{0},\mathbf{B}=\mathbb{1}_N$, the objective is simply absolute wealth.

- For $\mathbf{a}=\mathbf{0}, \mathbf{B}=\mathbb{1}_N - \frac{\mathbf{p}_T\boldsymbol{\beta}^\prime}{\boldsymbol{\beta}^\prime\mathbf{p}_T}$, for a secondary allocation $\boldsymbol{\beta}$, this represents an objective gauging that an allocation with outperform some reference allocation.

- For $\mathbf{a}=-\mathbf{p}_T, \mathbf{B}=\mathbb{1}_N$, we obtain net profits.

With prices being represented by some function of the market invariants $g$, we can get a glorified Taylor expansion called the *gamma* approximation to approximate the distribution of the objective:

$$ \Psi_{\boldsymbol{\alpha}}\approx\Theta_{\boldsymbol{\alpha}} = \theta_{\boldsymbol{\alpha}} + \boldsymbol{\Delta}^\prime_{\boldsymbol{\alpha}}\mathbf{X} + \frac{1}{2}\mathbf{X}^\prime\boldsymbol{\Gamma}_{\boldsymbol{\alpha}}\mathbf{X}$$

$$ \theta_{\boldsymbol{\alpha}} = \sum_{n=1}^N\alpha_n a_n + \sum_{n,m=1}^N\alpha_nB_{nm}g^{(m)}(\mathbf{0}) $$
$$ \boldsymbol{\Delta}_{\boldsymbol{\alpha}} = \sum_{n,m=1}^NB_{nm}\frac{\partial g^{(m)}}{\partial\mathbf{x}}|_{\mathbf{x}=\mathbf{0}}$$

$$\boldsymbol{\Gamma}_{\boldsymbol{\alpha}} = \sum_{n,m=1}^N\alpha_nB_{nm}\frac{\partial^2 g^{(m)}}{\partial\mathbf{x}\partial\mathbf{x}^\prime}|_{\mathbf{x}=\mathbf{0}} $$

If the invariants are normally distributed $\mathbf{X}\sim N(\boldsymbol{\mu},\boldsymbol{\Sigma})$, which the first part of this software can tell you, then the characteristic function of the approximate objective has a closed form expression in terms of all of these parameters, which makes computing the non-central moments of the approximate objective easy!! For example, $E[\Theta_{\boldsymbol{\alpha}}] = i^{-1}\phi^\prime_{\Theta_{\boldsymbol{\alpha}}}(0)$. These expressions are extremely long, and tedious, and won't be typed out by me right now!

So, this now gets us from estimating the distribution of market invariants, with a few different ways, and under some assumptions now being able to determine the distribution of an objective that takes a specific allocation of these invariants. Neat!

### Indices of satisfaction

I am not going to go into the entire theory of this, just the highlights. To say one allocation is "better" than another technically requires what's called strong or weak stochastic dominance. This is a huge field of itself, and so is the compromises we can come up with to instead evaluate allocations, such as indices of satisfaction. I'll cover a few here.

An investor will enjoy a generic outcome of an allocation $\boldsymbol{\alpha}$ by means of their utility function $u(\psi)$, assuming the realization of the outomce $\Psi_{\boldsymbol{\alpha}}=\psi$ actually occurs. It might then make sense to weight each outcome and get an expected utility from an allocation:

$$E[u(\Psi_{\boldsymbol{\alpha}})] = \sum_\mathbb{R}u(\psi)f_{\Psi_{\boldsymbol{\alpha}}}(\psi)d\psi$$

with $f$ the pdf of the objective (see above). This is the *von Neumann-Morgenstern* specification, but we want this to be in units of \$\$ so we instead consider the *certainty-equivalent*. This states the risk-free amount of moolah that would make someone as satisfied as the risky allocation:

$$CE(\boldsymbol{\alpha}) = u^{-1}(E[u(\Psi_{\boldsymbol{\alpha}})])$$

Because investors pursue the largest amount from their objectives, utility functions must be monotonically increasing, which implies that this inverse is always defined (and also increasing). If you wanted to get crazy with this, you can show that this utility function if the cumulative distribution function of the investors subjective before-the-fact hunch on the result of his investments :exploding_head:. That also reduces the certainty-equivalent to the quantile of this subjective distribution with a confidence level being the expected subjective grade.

