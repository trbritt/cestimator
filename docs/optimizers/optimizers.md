## Optimizer 

For the time being, we implement two types of portfolio optimization schemes: random allocations, and allocations
that lie in the "efficient frontier". For a given set of returns to consider $\mathbf{p_t}$, with mean $\bar{p}$, dispersed according to $\boldsymbol{\Sigma}$, the portfolio is uniquely described as a linear combination of assets $\mathbf{p_t} = \{p^{(0)}_t, p^{(1)}_t,\cdots, p^{(N)}_t\}$ via the coefficients of their linear combination $\mathbf{w}\in\mathcal{M}_{N\times 1}(\mathbb{R})$. In this case, the optimization, discovered by Markowitz, satisfies the following quadratic programming problem:

$$ 
\begin{cases}
\text{minimize} & -\bar{\mathbf{p}}^\prime \mathbf{w} + \mu \mathbf{w}^\prime \boldsymbol{\Sigma}\mathbf{w}&\\
\text{subject to} & \mathbf{1}_{1\times N}\mathbf{w} = 1, & \mathbf{w}\succeq 0
\end{cases} 
$$

In this case, we consider that all positions must be long, and that the sum of allocations must satisfy the budget equality (here normalized to 1). Many extensions are possible to implement, such as allowing short positions via the introduction of variables and constraints:

$$\mathbf{w}_\mathrm{long}\succeq 0,\quad\mathbf{w}_\mathrm{short}\succeq 0,\quad \mathbf{w}\equiv\mathbf{w}_\mathrm{long}-\mathbf{w}_\mathrm{short},\quad \mathbb{1}_{1\times N}\mathbf{w}_\mathrm{short}\leq \eta\mathbf{1}_{1\times N}\mathbf{w}_\mathrm{long} $$

where $\eta$ is the limiter for the short position at the beginning of the period to some fraction of the total long position at the beginning of the period. 

You can also easily include transaction costs. Assuming the buying and selling fees are $f_\mathrm{buy},f_\mathrm{sell}\geq 0$, we can make the change of variables using $\mathbf{u}_\mathrm{buy},\mathbf{u}_\mathrm{sell}$ which determine the amount of each asset we buy and sell before the holding period:
$$ \mathbf{w}\to\mathbf{w}+\mathbf{u}_\mathrm{buy}-\mathbf{u}_\mathbf{sell}, \quad \mathbf{u}_\mathrm{sell}\succeq 0,\quad \mathbf{u}_\mathrm{sell}\succeq 0$$

and the budget constraint now becomes the requirement taht the initial buying and selling nets zero cash:

$$(1-f_\mathrm{sell})\mathbf{1}_{1\times N}\mathbf{u}_\mathrm{sell} = (1+f_\mathrm{buy})\mathbb{1}_{1\times N}\mathbf{u}_\mathrm{buy} $$

To implement the solution to the quadratic problem, we use `eiquadprog`, originally created by Gabriele Buondonno, which can be found [here](https://www.cs.cmu.edu/~bstephe1/eiquadprog.hpp).

An example usage is given in `markowitz.cpp` in the `examples` directory.