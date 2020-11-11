Painleve I
========================================

Following the paper by Costin and Dunne [1904.11593](https://arxiv.org/abs/1904.11593).

ToDo
========================================
- Create the Series for the Painleve I
- Borel Transform


Creating the Series
----------------------------------------
$$y''(x) - 6 y^2(x) + x = 0$$

- want smooth real solution for large x
- Ecalle time $t = \frac{(24 x)^{5/4}}{30}$
- $y(x) = - \sqrt{\frac{x}{6}} ( 1 + h(t) )$
- now for $h(t)$ we have

$$\ddot{h}(t) + \frac{1}{t} \dot{h}(t) + h(t) \left( 1 + \frac{1}{2} h(t) \right)
- \frac{4}{25 t^2} \left(1 + h(t)\right) = 0$$

- we will expand
$$h(t) = \sum_{n=1}^\infty \frac{a_n}{t^{2n}}$$
- expand series to lowest needed order such that we have all terms to solve for the lowest order
- solve for lowest order
- expand again for next order
- not sure what it's good for, but I found that one can also write it analytically as
$$\sum_{n=1}^\infty \left[ \left( 2a_n ( 2n(n+1) -1 ) + a_{n+1} \right) \frac{1}{t^{2n+1}} + \frac{1}{2} \sum_{m=1}^n \frac{a_n a_{m-n+1}}{t^{2m+2}} \right] - \frac{4}{25t^2} + \frac{a_1}{t^2} = 0$$

From this one can immediately see $a_1=4/25$.

The algorithm is as follows:

- Start with an empty array `aArr` that will contain all known coefficients
- Loop over `order` starting from 1
	- Expand the DEQ to the needed order (so that everything to solve for `a[order]` is there). This might need some work with pencil and paper to figure out to which order one has to go.
	- Apply all coefficients of `aArr`.
	- Solve for `a[order]` and add it to `aArr`.

Possible improvements could be:

- maybe better to have the expansion step directly use the known coefficients and not applying them later.

Richardson's Extrapolation
========================================
We want to get the asymptotic form of $a_n$, $n\to\infty$.
I think the idea is to guess
$$a_n \sim C (-1)^{n+1} \Gamma\left( 2n - \frac{1}{2} \right)\,,$$
and then define new coefficients
$$\alpha_n \sim \frac{a_n}{(-1)^{n+1} \Gamma\left( 2n - \frac{1}{2} \right)}\,,$$
which should converge to $C$ which we want to find.
We will find it to higher order using Richardson's extrapoltion.

The theory is described in Bender, Orszag on p. 375.
One assumes that the error we make when evaluating a series $A_\infty=\sum_{k=0}^\infty a_k$ only to order $n$ goes like
$$ A_n = \sum_{k=0}^n \alpha_k \sim A_\infty + Q_1 n^{-1} + Q_2 n^{-2} + \ldots\,,\quad n\to \infty\,,$$
where we want to extract $A_\infty$ from only finitely many terms.
Using different lengths $n$ one can extract some of the $Q_i$ and thus improve the order.
For example if we use $A_n$ and $A_{n+1}$ we can extract $Q_1$ and thus improve convergence by
subtracting $Q_1/n$ from our result.

The general formula if we have $n+N$ terms and want to do a Richardson's extrapolation of
order $N$ is
$$A_\infty = \sum_{k=0}^N \frac{A_{n+k} (n+k)^N (-1)^{k+n}}{k!(N-k)!}$$
In my code I called the length of the series $n$ and thus I replaced $n\to n-N$ everywhere.

Tritronquee Solution
========================================

I don't know how to get this solution. The paper states it is the one with the tritronquee
initial conditions in some paper, but I could not find the initial conditions in this paper.
I now used the expansion given in Dunnes paper,
$$- 0.1875543083 - 0.3049055603 x + 0.1055298557 x^2 - 0.05229396374 x^3$$
but my solutions don't really match it for intermediate x.

Pade-Borel Transform
========================================
Borel and Pade
----------------------------------------

We transform the function $h(t)$
$$\mathcal{B}[h](p) = \sum_{n=1}^\infty \frac{a_n}{(2n-1)!} p^{2n-1}\,.$$
I think the $(2n-1)!$ is needed because of the $\Gamma(2n-1/2)$ in the growth of the coefficients.

Denoting the at $N$ truncated Borel transform as $\mathcal{B}_N[h](p)$ and its
\textit{Pade-Borel transform} is of the form
$$\mathcal{PB}_N[h](p) = \frac{P_{N-1}(p)}{Q_N(p)}\,.$$
We can find the poles of the Borel tranform by calculating the zeros of the $N$th order
polynomial $Q_N$.
Apparently the zeros of $P_{N-1}$ are also important for finding what kind of singularity we have.

Singularities
----------------------------------------
In our case one can find from Darboux's theorem (ToDo!) that the poles at
$p = \pm i$ are square root branch points and have the form
$$\mathcal{B}[h](p) \sim \frac{c}{\sqrt{p \mp i}}\,,\quad p\to \pm i\,.$$
