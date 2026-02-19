---
title: High-Dimensional Variance
date: 2025-12-28
layout: blog
cover: /assets/images/2025-12-28/whitening.png
tags:
  - Linear Algebra
intro: A useful view of a covariance matrix that is a natural generalization of variance to higher dimensions.
---

When $X$ is a random scalar, we call $\mathbb{V}[X]$ the _variance_ of the $X$, and define it as the average square distance from the mean:

$$
\begin{equation}
    \mathbb{V}[X]:=\sigma^2=\mathbb{E}\left[(X-\mathbb{E}[X])^2\right]
\end{equation}
$$

The mean $\mathbb{E}[X]$ is a measure of the _central tendency_ of $X$, while the variance $\mathbb{V}[X]$ is a measure of the _dispersion_ about this mean. Both metrics are moments of the distribution.

However, when $\textbf{X}$ is a _random vector_ or an ordered list of random scalars

$$
\begin{equation}
X=\begin{bmatrix}
X_1\\
X_2\\
\vdots\\
X_n
\end{bmatrix},
\end{equation}
$$

then, at least in my experience, people do not talk about the variance of $\textbf{X}$; rather, they talk about the _covariance matrix_ of $\textbf{X}$, denoted $\text{cov}(\textbf{X})$

$$
\begin{equation}
    \text{cov}(\mathbf{X}):=\Sigma=\mathbb{E}\left[(\mathbf{X}-\mathbb{E}[\mathbf{X}])(\mathbf{X}-\mathbb{E}[\mathbf{X}])^\top\right].
\end{equation}
$$

Now Equations $1$ and $3$ look quite similar, and I think it is natural to eventually ask: is a covariance matrix just a _high_- or _multi-dimensional variance_ of $\mathbf{X}$? Does this notation make sense?

$$
\begin{equation}
    \mathbb{V}[\mathbf{X}]:=\text{cov}[\mathbb{X}].
\end{equation}
$$

Of course, this is not my idea, but this was not hot I was taught to think about covariance matrices. But I think it is an interesting connection. So let's explore this idea that a covariance matrix is just high-dimensional variance.

#### Definitions

To start, it is clear from definitions that a covariance matrix is related variance. We can write the outer product in Equation $3$ explicitly as the following:

$$
\begin{equation}
X=\begin{bmatrix}
    \mathbb{E}[(\mathbf{X}_1-\mathbb{E}[\mathbf{X}_1])(\mathbf{X}_1-\mathbb{E}[\mathbf{X}_1])]&\ldots &\mathbb{E}[(\mathbf{X}_n-\mathbb{E}[\mathbf{X}_n])(\mathbf{X}_n-\mathbb{E}[\mathbf{X}_n])]\\
    \vdots &\ddots &\vdots\\
    \mathbb{E}[(\mathbf{X}_n-\mathbb{E}[\mathbf{X}_n])(\mathbf{X}_1-\mathbb{E}[\mathbf{X}_1])]&\ldots& \mathbb{E}[(\mathbf{X}_n-\mathbb{E}[\mathbf{X}_n])(\mathbf{X}_n-\mathbb{E}[\mathbf{X}_n])]\\
\end{bmatrix},
\end{equation}
$$

Clearly the diagonal element in $\Sigma$ are the variances of the scalars $X_i$. So
$\text{cov}[\mathbb{X}]$ still capture the dispersion of each $X_i$.

But what are the cross-terms? These are the _covariances_ of the pairwise combinations of $X_i$ and $X_j$, defined as

$$
\begin{equation}
    \sigma_{ij}:=\text{cov}\left(X_i,X_j\right)=\mathbb{E}\left[(X_i-\mathbb{E}[X_i])(X_j-\mathbb{E}[X_j])\right].
\end{equation}
$$

An obvious observation is that Equation $6$ is a genealisation of Equation $1$. Thus, variance $(i=j)$ is the simply a special case of covariance:

$$
\begin{equation}
    \mathbb{E}[X_i]=\text{cov}(X_i,X_i).
\end{equation}
$$

Furthermore, we can see that $\sigma_{ij}$ is not just a function of univariate variance of each random variable; it is also a function of whether the two variances are correlated with the other. As a simple thought experiment, imagine that $X_i$ and $X_j$ both had high variances but were uncorrelated with each other, meaning that there was no relationship between large (small) values of $X_i$ and large (small) value of $X_j$. Then we would still expect $\sigma_{ij}$ to be small (Figure $1$).

{% include widgets/figure.html num=1 src="/assets/images/2025-12-28/covariance.png" alt="Covariance scatter plots" caption="Empirical distributions of two random variables $X_1$ and $X_2$. Each subplot is labeled with three numbers representing $\sigma_1,\sigma_2$, and $\sigma_{12}$ respectively. (Black subplots) $X_1$ and $X_2$ are negatively correlated with Pearson correlation coefficient $\rho=-0.7$. (Blue subplots) $X_1$ and $X_2$ are uncorrelated. (Green subplots) $X_1$ and $X_2$ are positively correlated $(\rho=0.7)$." %}

This is the intuitive reason that $\sigma_{ij}$ can be written in terms of the \*Pearson correlation coefficient $\rho_{ij}$

$$
\begin{equation}
    \sigma_{ij}=\rho_{ij}\sigma_{i}\sigma_{j}
\end{equation}
$$

Again, the generalises Equation $1$ but with $\rho_{ij}=1$.

Using the $\sigma_{ij}$ notation from Equation $8$ we can rewrite the covariance matrix as

$$
\begin{equation}
    \Sigma=
    \begin{bmatrix}
        \sigma_1^2 & \rho_{12}\sigma_1\sigma_2 & \ldots & \rho_{1n}\sigma_1\sigma_n \\
        \rho_{21}\sigma_2\sigma_2 & \sigma_2^2 & \ldots & \rho_{2n}\sigma_2\sigma_n \\
        \vdots & \ddots & \ldots &  \vdots \\
        \rho_{n1}\sigma_n\sigma_1 & \rho_{n2}\sigma_2\sigma_n & \ldots & \sigma_n^2\\
    \end{bmatrix}.
\end{equation}
$$

And with a little algebra, we can decompose the matrix in Equation $9$ into a form that looks very similar to a multi-dimensional Equation $8$.

$$
\begin{equation}
    \Sigma=
    \begin{bmatrix}
        \sigma_1 &&&\\
        &\sigma_2 &&\\
        &&\ddots& \\
        &&& \sigma_n \\
    \end{bmatrix}
    \begin{bmatrix}
        1 & \rho_{12} & \ldots & \rho_{1n}\\
        \rho_{21} & 1 & \ldots & \rho_{2n}\\
        \ldots & \ldots & \ddots & \ldots\\
        \rho_{n1} & \rho_{n2} & \ldots & 1\\
    \end{bmatrix}
    \begin{bmatrix}
        \sigma_1 &&&\\
        &\sigma_2 &&\\
        &&\ddots& \\
        &&& \sigma_n \\
    \end{bmatrix}.
\end{equation}
$$

The middle matrix is _correlation matrix_, which captures the Pearson (linear) correlation between all the variables in $\mathbf{X}$. So a covariance matrix of $\mathbf{X}$ captures both the dispersion of each $X_i$ and how it covaries with other random variables. If we think of scalar variance as a univariate covariance matrix, then one-dimensional case can be written as

$$
\begin{equation}
    \Sigma=[\text{cov}(X,X)]=[\sigma][1][\sigma].
\end{equation}
$$

Equation $11$ is not useful from a practical point of view, but in my mind, it underscores the view that univariate variance is a special case of covariance (matrices).

Finally, note that Equation $10$ gives us a useful way to compute the correlation matrix $\mathbf{C}$ from the covariance matrix $\Sigma$. It is

$$
\begin{equation}
    \mathbf{C}:=
    \begin{bmatrix}
        1/\sigma_1 & &  \\
                   &\ddots& \\
                   & & 1/\sigma_n \\
    \end{bmatrix}
    \Sigma
    \begin{bmatrix}
        1/\sigma_1 & &  \\
                   &\ddots& \\
                   & & 1/\sigma_n \\
    \end{bmatrix}
\end{equation}
$$

This is because the inverse of a diagonal matrix is simply the reciprocal of each element along the diagonal.

#### Properties

Some important properties of univariate variance have multidimensional analogs. For example, let $a$ be a non-random number. Recall that the variance of $aX$ scales with $a^2$ or that

$$
\begin{equation}
    \mathbb{V}[aX]=a^2\mathbb{V}[X]=a^2\sigma^2.
\end{equation}
$$

And in the general case, $a$ is a non-random matrix $\mathbf{A}$ and

$$
\begin{equation}
\begin{aligned}
    \mathbb{V}[\mathbf{AX}]&=\mathbf{A}^\top\Sigma\mathbf{A} \\
                         &=\mathbb{E}\left[(\mathbf{AX}-\mathbb{E}[\mathbf{AX}])(\mathbf{AX}-\mathbb{E}[\mathbf{AX}])^\top\right]\\
                         &=\mathbf{A}\mathbb{E}\left[(\mathbf{X}-\mathbb{E}[\mathbf{X}])(\mathbf{X}-\mathbb{E}[\mathbf{X}])^\top\right]\mathbf{A}^\top \\
                         &=\mathbf{A}\Sigma\mathbf{A}^\top
\end{aligned}
\end{equation}
$$

So both $\mathbb{V}[\mathbf{X}]$ and $\mathbb{V}[\mathbf{X}]$ are quadratic with respect to multiplicative constants!

As a second example, recall that the univariate variance of $a+X$ is just variance of $X$:

$$
\begin{equation}
    \mathbb{V}[a+X]=\mathbb{V}[X].
\end{equation}
$$

This is intuitive. A constant shift in the distribution does not change its dispersion. And again, in the general case, a constant shift in a multivariate distribution does not change its dispersion:

$$
\begin{equation}
\begin{aligned}
    \mathbb{V}[\mathbf{A}+\mathbf{X}]&=\mathbb{E}\left[(\mathbf{A}+\mathbf{X}-\mathbb{E}[\mathbf{A}+\mathbf{X}])(\mathbf{A}+\mathbf{X}-\mathbb{E}[\mathbf{A}+\mathbf{X}])^\top\right]\\
    &=\mathbb{E}\left[(\mathbf{X}-\mathbb{E}[\mathbf{X}])(\mathbf{X}-\mathbb{E}[\mathbf{X}])^\top\right]\\
    &=\Sigma
\end{aligned}
\end{equation}
$$

Finally, a standard decomposition of variance is to write it as

$$
\begin{equation}
    \mathbb{E}[X]=\mathbb{E}\left[X^2\right]-\mathbb{E}[X]^2
\end{equation}
$$

And this standard decomposition has a multi-dimensional analog:

$$
\begin{equation}
\begin{aligned}
    \mathbb{V}[\mathbf{X}]&=\mathbb{E}\left[(\mathbf{X}-\mathbb{E}[\mathbf{X}])(\mathbf{X}-\mathbb{E}[\mathbf{X}])^\top\right]\\
    &=\mathbb{E}\left[\mathbf{X}\mathbf{X}^\top-2\mathbb{E}[\mathbf{X}]\mathbf{X}^\top+\mathbb{E}[\mathbf{X}]\mathbb{E}[\mathbf{X}]^\top \right]\\
    &=\mathbb{E}\left[\mathbf{X}\mathbf{X}^\top\right]-\mathbb{E}[\mathbf{X}]\mathbb{E}[\mathbf{X}]^\top
\end{aligned}
\end{equation}
$$

Regardless, my point isnot to provide detailed or comprehensive proofs, but only to underscore that covariance matrices have properties that indicate they are simply high-dimensional (co)variances.

#### Non-negativity

Another neat connection is that covariance matrices are **positive semi-definite** (PSD), which I'll
denote with

$$
\begin{equation}
    \mathbb{V}[X]=\Sigma\geq 0,
\end{equation}
$$

while univariate variance are non-negative number:

$$
\begin{equation}
    \mathbb{V}[X]=\sigma^2\geq 0.
\end{equation}
$$

In this view, the _Cholesky decomposition_ of a PSD matrix is simply a high-dimensional square root! So the Cholesky factor $\mathbf{L}$ can be viewed as the high-dimensional standard deviation of $X$, since

$$
\begin{equation}
    \Sigma=\mathbf{L}\mathbf{L}^\top\implies\mathbf{L}\approx\sigma
\end{equation}
$$

Here, I am abusing notation to use $\approx$ to mean "analogous to".

#### Precision and whitening

The precision of a scalar random variable is the reciprocal of its variance:

$$
\begin{equation}
    p=\frac{1}{\sigma^2}
\end{equation}
$$

Hopefully the name is somewhat intuitive. When a random variable has high precision, this means it has low variance and thus a smaller range of possible outcomes. A common place that precision arises is when _whitening_ data, or standardising it to have zero mean and unit variance. We can do this by subtracting the mean $X$ and then dividing it by it standard deviation:

$$
\begin{equation}
    Z=\frac{X-\mathbb{E}[X]}{\sigma}
\end{equation}
$$

This is sometimes referred to as _z-scoring_.

What's the multivariate analog to this? We can define the precision matrix $\mathbf{P}$ as the inverse of the covariance matrix or

$$
\begin{equation}
    \mathbf{P}=\Sigma^{-1}.
\end{equation}
$$

Given the Colesky decomposition in Equation $21$ above, we can compute the Cholesky factor of the precision matrix as:

$$
\begin{equation}
\begin{aligned}
    \mathbf{P}&=\left(\mathbf{L}\mathbf{L}^{\top}\right)^{-1}\\
              &=\left(\mathbf{L}^{\top}\right)^{-1}\mathbf{L}^{-1}\\
              &=\left(\mathbf{L}^{-1}\right)^{\top}\mathbf{L}^{-1}\\
\end{aligned}
\end{equation}
$$

So the multivariate analog to Equation $$ is

$$
\begin{equation}
    \mathbf{Z}=\mathbf{L}^{-1}(\mathbf{X}-\mathbb{E}[\mathbf{X}]).
\end{equation}
$$

The geometric or visual effect of this operation is to apply a linear transformation (rotation) of our data (samples of the random vector) with covariance matrix $\Sigma$ into a new set of variables with an identity covariance matrix (Figure $2$).

{% include widgets/figure.html num=2 src="/assets/images/2025-12-28/whitening.png" alt="whitening" caption="(Left) Samples of a random variable $X$ with variances $\sigma_1=1.5$ and $\sigma_2=0.5$ and correlation $\rho=0.5$. (Right) Samples rotated by the Cholesky factor of the known covariance matrix, as in Equation $26$. The effect is that the rotated samples are uncorrelated with unit variance." %}

Ignoring the mean, we can easily verify this transformation works:

$$
\begin{equation}
\begin{aligned}
\mathbb{V}[Z]&=\mathbb{E}\left[\left(\mathbf{L}^{-1}\mathbf{X}\right)\left(\mathbf{L}^{-1}\mathbf{X}\right)^\top\right]\\
        &=\mathbb{E}\left[\mathbf{L}^\top\mathbf{X}\mathbf{X}^\top\left(\mathbf{L}^{-1}\right)^\top\right]\\
        &=\mathbf{L}^{-1}\mathbb{E}\left[\mathbf{X}\mathbf{X}^\top\right]\left(\mathbf{L}^{-1}\right)^\top\\
        &=\mathbf{L}^{-1}\Sigma\left(\mathbf{L}^\top\right)^{-1}\\
        &=\mathbf{I}\\
\end{aligned}
\end{equation}
$$

This derivation is a little more tedious with a mean, but hopefully the idea is clear. _Why_ the Cholesky decomposition actually works is a deeper idea, one worth its own post, but I think the essential ideas are also captured in **principal component analysis** (PCA).

#### Summary statistic

It would be useful to summarise the information in a covariance matrix with a single number. To my knowledge, there are at least two such summary statistic that capture different types of information.

**Total variance.** The _total variance_ of a random vector $\mathbf{X}$ is the trace of its covariance matrix or

$$
\begin{equation}
    \sigma_{\text{tv}}^2:=\text{tr}(\Sigma)=\sum_{i=1}^n\sigma_{i}^2.
\end{equation}
$$

We can see that variance is a scalar that summarises the variance across the components of $\mathbf{X}$. The concept is used in PCA, where total variance is preserved across the transformation. Of course, in the-dimensional case, total variance is simply variance, $\sigma_{\text{tv}}^2=\sigma^2$.

**Generlised variance.** The generalised variance ([Wilks, 1932](URL)) of a random vector $\mathbf{X}$ is the determinant of its covariance matrix or

$$
\begin{equation}
    \sigma_{\text{gv}}^2:=\text{det}(\Sigma)=|\Sigma|.
\end{equation}
$$

There are nice [geometric](https://www.youtube.com/watch?v=Ip3X9LOh2dk) [interpretations]() of the determinant, but prehaps the simplest way to think about it here is that the determinant is equal to the product of the eigenvalues of $\Sigma$ or

$$
\begin{equation}
    |\Sigma|=\prod_{i=1}^n\lambda_i.
\end{equation}
$$

So we can think of generalised variance as capturing the magnitude of the linear
transformation represented by $\Sigma$.

Generalised variance is quite different from total variance. For example, consider the
two-dimensional covariance matrix implied by the following values:

$$
\begin{equation}
    \sigma_1=\sigma_2=\rho=1
\end{equation}
$$

Clearly, the total variance is two, while the determinant is one. Now imagine that the variables are highly correlated, that $\rho=0.98$. Then the total variance is still $2$, but the matrix becomes "more singular" (it is roughly $0.039$). So total variance, as its name suggests, really just summarises the dipersion in $\Sigma$, while generalised variance also captures how the variables in $\mathbf{X}$ covary. When the variable in $\mathbf{X}$ are highly (un-) correlated, generalised variance will be low (high).

#### Example

Let's end with two interesting examples that use the ideas in this post.

**Multivariate normal.** First, recall that the probabiliy density function (PDF) for a standard normal random variable is

$$
\begin{equation}
    p\left(x;\mu,\sigma^2\right)=\frac{1}{\sqrt{2\pi}\sigma}\exp\left\lbrace -\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2 \right\rbrace.
\end{equation}
$$

We can immediately see that the squared term is just the square of Equation $23$. Now armed the interpretation that a covariance matrix is high-dimensional variance, consider the PDF for a multivariate normal random variable:

$$
\begin{equation}
    p\left(\mathbf{x};\mu,\Sigma \right)=\frac{1}{(2\pi)^{-n/2}|\Sigma|^{-1/2}}\exp\left\lbrace -\frac{1}{2}(\mathbf{x}-\mu)^{\top}\Sigma^{-1}(\mathbf{x}-\mu) \right\rbrace.
\end{equation}
$$

We can see that the [Mahalanobis distance](https://en.wikipedia.org/wiki/Mahalanobis_distance) is a multivariate whitening. And the variance is the normalising term,

$$
\begin{equation}
    \frac{1}{\sqrt{2\pi}\sigma}
\end{equation}
$$

has a multivariate analog that is generalised variance above.

**Correlated random variables.** Consider a scalar random variable $Z$ with unit variance, $\mathbb{V}[Z]=1$. We can transform this into a random variable $X$ with variance $\sigma^2$ by multiplying $Z$ by $\sigma$:

$$
\begin{equation}
    \mathbb{V}[\sigma Z]=\sigma^2\mathbb{V}[Z]=\sigma^2.
\end{equation}
$$

What is the mutli-dimensional version of this? If we have two random variables

$$
\begin{equation}
    Z=
    \begin{bmatrix}
        Z_1 \\
        Z_2 \\
    \end{bmatrix},
\end{equation}
$$

how can we transform them into a random variable $\mathbf{X}$ with covariance matrix $\Sigma$? Clearly, we multiply $\mathbf{Z}$ by the Cholesky factor, a la Equation $26$.

What this suggest, however, is a generic algorithm for generating correlated random variable: we multiply $\mathbf{Z}$ by the Cholesky factor of the covariance matrix when $\sigma_1=\ldots=\sigma_n=1$. In the two-by-two case, that's

$$
\begin{equation}
    \text{cholesky}\left(\begin{bmatrix}1&\rho\\ \rho&1\\\end{bmatrix} \right)=\mathbf{L}\mathbf{L}^{\top}=\begin{bmatrix}1&0\\ \rho&\sqrt{1-\rho^2}\end{bmatrix} \begin{bmatrix}1&\rho\\ 0&\sqrt{1-\rho^2} \end{bmatrix}
\end{equation}
$$

This suggest that an algorithm: draw two i.i.d random variables $Z_1$ and $Z_2$ both with unit variance. The set

$$
\begin{equation}
\begin{aligned}
    X_1&:=Z_1\\
    X_1&:=Z_1\rho+Z_2\sqrt{1-\rho^2}\\
\end{aligned}
\end{equation}
$$

Of course, this is nice because it can be vectorized and extends to an arbitrary number of random variables. And we can easily account for non-unit variances if we could like.

#### Conclusion

With the proper framing, it is fairly natural to think of $\Sigma$ as simply variance and to denote it as $\mathbb{V}[\mathbf{X}$. This framing is useful because makes certain properties of covariance matrices almost obvious, such as why they are positive semi-definiteness or why their inverses appear when whitening data. However, high-dimensional variance has properities that are not important in one-dimension, such as the correlation between the variables in $\mathbf{X}$. Thus in my mind, the best framing is that univeriate variance is really a special case of covariance matrices. However, in my mind, either reframing is useful for gaining a deeper intuition for the material.
