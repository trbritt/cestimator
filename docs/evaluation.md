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

If the invariants are normally distributed $\mathbf{X}\sim N(\boldsymbol{\mu},\boldsymbol{\Sigma})$, which the first part of this software can tell you, then the characteristic function of the approximate objective has a closed form expression in terms of all of these parameters, which makes computing the non-central moments of the approximate objective easy!! For example, $E[\Theta_{\boldsymbol{\alpha}}] = i^{-1}\phi^\prime_{\Theta_{\boldsymbol{\alpha}}}(0)$. 

#### The characteristic function

Defining $\mathbf{V}=\boldsymbol{\Gamma}\boldsymbol{\Sigma}$ and:

$$b_{\boldsymbol{\alpha}} = \theta_{\boldsymbol{\alpha}}+\boldsymbol{\Delta}_{\boldsymbol{\alpha}}^\prime \boldsymbol{\mu} + \frac{1}{2}\boldsymbol{\mu}^\prime\boldsymbol{\Gamma}_{\boldsymbol{\alpha}}\boldsymbol{\mu}$$

$$\mathbf{w}_{\boldsymbol{\alpha}}=\boldsymbol{\Delta}_{\boldsymbol{\alpha}}+\boldsymbol{\Gamma}_{\boldsymbol{\alpha}}\boldsymbol{\mu}$$


you can show, through a massive headache and algebra, that the characteristic function of the approximate objective distribution is:

$$ \phi_{\Theta_{\boldsymbol{\alpha}}}(\omega) = v^{-1/2}e^u$$

$$ u(\omega) = i\omega b - \frac{1}{2}\mathbf{w}^\prime\boldsymbol{\Sigma}(\mathbb{1}-i\omega\mathbf{V})^{-1}\mathbf{w}$$

$$ v(\omega) = |\mathbb{1}-i\omega\mathbf{V}|$$

This now lets you take the derivatives of the characteristic function in terms of u (enter Wolfram Alpha), and allows you to put it in terms of everything you already know, albeit with some headache of symbol vomit on your screen. 

Punchline: you estimated the distribution of invariants, aka obtained $\boldsymbol{\mu}$ and $\boldsymbol{\Sigma}$. You can now therefore choose a combination of these invariants, and choose the benchmark (objective) against which you will rate different allocations. This will allow you determine directly the characteristic function of your objective (equivalent to PDF and CDF), and therefore its non-central moments. These non-central moments are then the final ingredients needed to evaluate various metrics of quality of the objective for that allocation, which we show ...

---


### Indices of satisfaction

... right now! I am not going to go into the entire theory of this, just the highlights. To say one allocation is "better" than another technically requires what's called strong or weak stochastic dominance. This is a huge field of itself, and so is the compromises we can come up with to instead evaluate allocations, such as indices of satisfaction. I'll cover a few here.

An investor will enjoy a generic outcome of an allocation $\boldsymbol{\alpha}$ by means of their utility function $u(\psi)$, assuming the realization of the outomce $\Psi_{\boldsymbol{\alpha}}=\psi$ actually occurs. It might then make sense to weight each outcome and get an expected utility from an allocation:

$$E[u(\Psi_{\boldsymbol{\alpha}})] = \sum_\mathbb{R}u(\psi)f_{\Psi_{\boldsymbol{\alpha}}}(\psi)d\psi$$

with $f$ the pdf of the objective (see above). This is the *von Neumann-Morgenstern* specification, but we want this to be in units of \$\$ so we instead consider the *certainty-equivalent*. This states the risk-free amount of moolah that would make someone as satisfied as the risky allocation:

$$CE(\boldsymbol{\alpha}) = u^{-1}(E[u(\Psi_{\boldsymbol{\alpha}})])$$

Because investors pursue the largest amount from their objectives, utility functions must be monotonically increasing, which implies that this inverse is always defined (and also increasing). If you wanted to get crazy with this, you can show that this utility function if the cumulative distribution function of the investors subjective before-the-fact hunch on the result of his investments :exploding_head:. That also reduces the certainty-equivalent to the quantile of this subjective distribution with a confidence level being the expected subjective grade.

Instead of gauging the risk of a given indec of satisfaction globally (certainty-equivalent is risk prone iff $u$ is convex or neutral iff $u$ is linear, etc), you'll be more inclined to be cautious of pursuing new gains, but wouldn't want to cut your losses in hopes of recovery (aka prospect theory), so we need a local measure of risk, which is indeed the *Arrow-Pratt absolute risk aversion*:

$$ A(\psi) = -\frac{\mathcal{D}^2u(\psi)}{\mathcal{D}u(\psi)}$$

with $\mathcal{D}$ the derivative operator. It is worth mentioning that it is more practical to parametrize $A$ instead of $u$, which is typically done by a set $(\gamma, \zeta, \eta)$ via:
$$ A(\psi) = \frac{\psi}{\gamma\psi^2 + \zeta\psi + \eta}$$ 
that encompasses a massive range of use cases. This let's you make the approximation:

$$CE(\boldsymbol{\alpha}) \approx E[\Psi_{\boldsymbol{\alpha}}] - \frac{A(E[\Psi_{\boldsymbol{\alpha}}])}{2}Var[\Psi_{\boldsymbol{\alpha}}]$$

where chances are you don't know these moments exactly, but you can compute now by estimation from the gamma approximation above and the location and dispersion estimators you computed using the `Cestimator::Estimator` classes ;)

You could instead also likewise enforce the idea that your losses should not exceed some threshold $L$ with some level of confidence $c$:

$$\mathbb{P}[w_T-W_{T+\tau} < L] \ge c$$

which leads you to the value at risk index of satisfaction:

$$VaR_c(\boldsymbol{\alpha}) = - Q_{\Psi_{\boldsymbol{\alpha}}}(1-c)$$

For a generic allocation, the Cornish-Fisher expansion lets you evaluate the quantile in terms of the quantile of the standard normal distribution $z(p)=\sqrt{2}\mathrm{erf}^{-1}(2p-1)$:

$$ Q_{\Psi_{\boldsymbol{\alpha}}}(1-c) \approx A_{\boldsymbol{\alpha}} + B_{\boldsymbol{\alpha}}z(1-c)+C_{\boldsymbol{\alpha}}z^2(1-c)$$

$$A = E[\Psi_{\boldsymbol{\alpha}}] - \frac{E[\Psi_{\boldsymbol{\alpha}}^3]-3E[\Psi_{\boldsymbol{\alpha}}^2]E[\Psi_{\boldsymbol{\alpha}}]+2E[\Psi_{\boldsymbol{\alpha}}]^3}{6\left(E[\Psi_{\boldsymbol{\alpha}}^2] - E[\Psi_{\boldsymbol{\alpha}}]^2\right)}$$

$$B=\sqrt{E[\Psi_{\boldsymbol{\alpha}}^2]-E[\Psi_{\boldsymbol{\alpha}}]^2}$$

$$C = \frac{E[\Psi_{\boldsymbol{\alpha}}^3] - 3E[\Psi_{\boldsymbol{\alpha}}^2]E[\Psi_{\boldsymbol{\alpha}}]+2E[\Psi_{\boldsymbol{\alpha}}]^3}{6(E[\Psi_{\boldsymbol{\alpha}}^2]-E[\Psi_{\boldsymbol{\alpha}}]^2)}$$

with these terms being evaluated via derivatives of the characteristic function:

$$i^kE[\Psi_{\boldsymbol{\alpha}}^k] = \frac{d^k\phi_{\Psi_{\boldsymbol{\alpha}}}(\omega)}{d\omega^k}\bigg|_{\omega=0}$$