---
title: High-Dimensional Variance
date: 2025-12-28
layout: blog
cover: /assets/images/2025-12-28/whitening.png
tags:
  - Linear Algebra
  - Probability and Statistics
intro: A useful view of a covariance matrix that is a natural generalization of variance to higher dimensions.
---

When $X$ is a random scalar, we refer to $\mathbb{V}[X]$ as the _variance_ of $X$, defining it as the expected squared distance from the mean:

$$
\begin{equation}
    \mathbb{V}[X]:=\sigma^2=\mathbb{E}\left[(X-\mathbb{E}[X])^2\right]
\end{equation}
$$

The mean $\mathbb{E}[X]$ describes the _central tendency_ of $X$, whereas the variance $\mathbb{V}[X]$ quantifies the _dispersion_ around this mean. Both are moments of the distribution.

However, when $\mathbf{X}$ is a _random vector_—that is, an ordered list of random scalars,

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

people, at least in my experience, tend not to speak of the variance of $\mathbf{X}$ as such; instead, they refer to the _covariance matrix_ of $\mathbf{X}$, denoted $\text{cov}(\mathbf{X})$:

$$
\begin{equation}
    \text{cov}(\mathbf{X}):=\Sigma=\mathbb{E}\left[(\mathbf{X}-\mathbb{E}[\mathbf{X}])(\mathbf{X}-\mathbb{E}[\mathbf{X}])^\top\right].
\end{equation}
$$

Equations $1$ and $3$ now appear rather similar, leading naturally to the question: is the covariance matrix simply a _higher- or multi-dimensional variance_ of $\mathbf{X}$? Does writing

$$
\begin{equation}
    \mathbb{V}[\mathbf{X}]:=\text{cov}[\mathbf{X}]
\end{equation}
$$

make sense?

Of course, this perspective is not original to me; it was not the way I was initially taught to think about covariance matrices. Nevertheless, I find the connection insightful. Let us therefore explore the idea that a covariance matrix is simply a high-dimensional generalisation of variance.

#### Definitions

From the definitions, it is evident that the covariance matrix is intimately related to variance. We can write the outer product in Equation $3$ explicitly as follows:

$$
\begin{equation}
\Sigma=\begin{bmatrix}
    \mathbb{E}[(X_1-\mathbb{E}[X_1])^2] & \ldots & \mathbb{E}[(X_1-\mathbb{E}[X_1])(X_n-\mathbb{E}[X_n])]\\
    \vdots & \ddots & \vdots\\
    \mathbb{E}[(X_n-\mathbb{E}[X_n])(X_1-\mathbb{E}[X_1])] & \ldots& \mathbb{E}[(X_n-\mathbb{E}[X_n])^2]\\
\end{bmatrix},
\end{equation}
$$

Clearly, the diagonal elements of $\Sigma$ are the variances of the scalar components $X_i$, so $\text{cov}[\mathbf{X}]$ still reflects the dispersion of each $X_i$.

What, then, are the off-diagonal entries? These are the _covariances_ between pairs $X_i$ and $X_j$, defined as

$$
\begin{equation}
    \sigma_{ij}:=\text{cov}\left(X_i, X_j\right)=\mathbb{E}\left[(X_i-\mathbb{E}[X_i])(X_j-\mathbb{E}[X_j])\right].
\end{equation}
$$

It is clear that Equation $6$ generalises Equation $1$, with the variance case being the special instance $i = j$:

$$
\begin{equation}
    \mathbb{V}[X_i]=\text{cov}(X_i, X_i).
\end{equation}
$$

Moreover, $\sigma_{ij}$ is affected not only by the univariate variances of the two variables but also by how much those variables are correlated with one another. For example, suppose $X_i$ and $X_j$ both exhibit high variance, but are uncorrelated—that is, there is no statistical relationship between large (or small) values of $X_i$ and those of $X_j$. In that case, we would expect $\sigma_{ij}$ to be small (as shown in Figure $1$).

{% include widgets/figure.html num=1 src="/assets/images/2025-12-28/covariance.png" alt="Covariance scatter plots" caption="Empirical distributions of two random variables $X_1$ and $X_2$. Each subplot is labelled with three numbers representing $\sigma_1,\sigma_2$, and $\sigma_{12}$ respectively. (Black subplots) $X_1$ and $X_2$ are negatively correlated with Pearson correlation coefficient $\rho=-0.7$. (Blue subplots) $X_1$ and $X_2$ are uncorrelated. (Green subplots) $X_1$ and $X_2$ are positively correlated $(\rho=0.7)$." %}

This intuitive argument is why $\sigma_{ij}$ can be written in terms of the Pearson correlation coefficient, $\rho_{ij}$:

$$
\begin{equation}
    \sigma_{ij}=\rho_{ij}\sigma_{i}\sigma_{j}
\end{equation}
$$

This again generalises Equation $1$, with $\rho_{ij}=1$ in the variance case.

Using the $\sigma_{ij}$ notation from Equation $8$, we can rewrite the covariance matrix as

$$
\begin{equation}
    \Sigma=
    \begin{bmatrix}
        \sigma_1^2 & \rho_{12}\sigma_1\sigma_2 & \cdots & \rho_{1n}\sigma_1\sigma_n \\
        \rho_{21}\sigma_2\sigma_1 & \sigma_2^2 & \cdots & \rho_{2n}\sigma_2\sigma_n \\
        \vdots & \ddots & \cdots &  \vdots \\
        \rho_{n1}\sigma_n\sigma_1 & \rho_{n2}\sigma_n\sigma_2 & \cdots & \sigma_n^2\\
    \end{bmatrix}.
\end{equation}
$$

A little algebra allows us to decompose the matrix in Equation $9$ as follows, resembling a multi-dimensional generalisation:

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

The central matrix is the _correlation matrix_, encoding the Pearson (linear) correlations between each pair of variables in $\mathbf{X}$. Thus, the covariance matrix of $\mathbf{X}$ captures both the individual dispersions of the $X_i$ and their covariances with the other variables. If we regard scalar variance as a univariate covariance matrix, the one-dimensional case may be written as

$$
\begin{equation}
    \Sigma=[\text{cov}(X,X)]=[\sigma][1][\sigma].
\end{equation}
$$

While Equation $11$ is of little practical use, it emphasises the idea that univariate variance is simply a special case of covariance matrices.

Equation $10$ also provides a straightforward method for deriving the correlation matrix $\mathbf{C}$ from the covariance matrix $\Sigma$:

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

This holds because the inverse of a diagonal matrix consists of the reciprocal of each diagonal entry.

#### Properties

Many important properties of univariate variance have direct analogues in multiple dimensions. For instance, if $a$ is a non-random scalar, recall that the variance of $aX$ scales by $a^2$:

$$
\begin{equation}
    \mathbb{V}[aX]=a^2\mathbb{V}[X]=a^2\sigma^2.
\end{equation}
$$

In the general case, with a non-random matrix $\mathbf{A}$,

$$
\begin{equation}
\begin{aligned}
    \mathbb{V}[\mathbf{A}\mathbf{X}]&=\mathbf{A}\Sigma\mathbf{A}^\top \\
                         &=\mathbb{E}\left[(\mathbf{A}\mathbf{X}-\mathbb{E}[\mathbf{A}\mathbf{X}])(\mathbf{A}\mathbf{X}-\mathbb{E}[\mathbf{A}\mathbf{X}])^\top\right]\\
                         &=\mathbf{A}\mathbb{E}\left[(\mathbf{X}-\mathbb{E}[\mathbf{X}])(\mathbf{X}-\mathbb{E}[\mathbf{X}])^\top\right]\mathbf{A}^\top \\
                         &=\mathbf{A}\Sigma\mathbf{A}^\top
\end{aligned}
\end{equation}
$$

Thus, both $\mathbb{V}[\mathbf{X}]$ and its multivariate transformations are quadratic in multiplicative constants.

As another example, recall that the univariate variance of $a+X$ is simply the variance of $X$:

$$
\begin{equation}
    \mathbb{V}[a+X]=\mathbb{V}[X].
\end{equation}
$$

Clearly, shifting a distribution by a constant does not affect its dispersion. In the multivariate case, a constant shift similarly leaves the covariance unchanged:

$$
\begin{equation}
\begin{aligned}
    \mathbb{V}[\mathbf{A}+\mathbf{X}]&=\mathbb{E}\left[(\mathbf{A}+\mathbf{X}-\mathbb{E}[\mathbf{A}+\mathbf{X}])(\mathbf{A}+\mathbf{X}-\mathbb{E}[\mathbf{A}+\mathbf{X}])^\top\right]\\
    &=\mathbb{E}\left[(\mathbf{X}-\mathbb{E}[\mathbf{X}])(\mathbf{X}-\mathbb{E}[\mathbf{X}])^\top\right]\\
    &=\Sigma
\end{aligned}
\end{equation}
$$

Further, variance has a standard decomposition:

$$
\begin{equation}
    \mathbb{V}[X]=\mathbb{E}\left[X^2\right]-\mathbb{E}[X]^2
\end{equation}
$$

The analogous result for vectors is

$$
\begin{equation}
\begin{aligned}
    \mathbb{V}[\mathbf{X}]&=\mathbb{E}\left[(\mathbf{X}-\mathbb{E}[\mathbf{X}])(\mathbf{X}-\mathbb{E}[\mathbf{X}])^\top\right]\\
    &=\mathbb{E}\left[\mathbf{X}\mathbf{X}^\top-2\mathbb{E}[\mathbf{X}]\mathbf{X}^\top+\mathbb{E}[\mathbf{X}]\mathbb{E}[\mathbf{X}]^\top \right]\\
    &=\mathbb{E}\left[\mathbf{X}\mathbf{X}^\top\right]-\mathbb{E}[\mathbf{X}]\mathbb{E}[\mathbf{X}]^\top
\end{aligned}
\end{equation}
$$

My aim here is not to provide exhaustive proofs, but rather to highlight that the covariance matrix retains properties analogous to 'variance', but in higher dimensions.

#### Non-negativity

A further illuminating connection is that covariance matrices are **positive semi-definite** (PSD):

$$
\begin{equation}
    \mathbb{V}[X]=\Sigma\geq 0,
\end{equation}
$$

whilst univariate variances are non-negative numbers:

$$
\begin{equation}
    \mathbb{V}[X]=\sigma^2\geq 0.
\end{equation}
$$

Accordingly, the _Cholesky decomposition_ of a PSD matrix is a natural generalisation of the scalar square root. The Cholesky factor $\mathbf{L}$ thus plays the role of a multi-dimensional standard deviation, since

$$
\begin{equation}
    \Sigma=\mathbf{L}\mathbf{L}^\top\implies\mathbf{L}\approx\sigma
\end{equation}
$$

Here, ‘$\approx$’ is used to indicate analogy.

#### Precision and whitening

The precision of a scalar random variable is simply the inverse of its variance:

$$
\begin{equation}
    p=\frac{1}{\sigma^2}
\end{equation}
$$

The intuition is that higher precision corresponds to lower variance, i.e. less spread in possible outcomes. Precision naturally arises when ‘whitening’ data, which involves standardising it to have zero mean and unit variance. This is achieved by subtracting the mean and dividing by the standard deviation:

$$
\begin{equation}
    Z=\frac{X-\mathbb{E}[X]}{\sigma}
\end{equation}
$$

This is sometimes called _z-scoring_.

What is the multivariate counterpart? Here, the precision matrix $\mathbf{P}$ is defined as the inverse of the covariance matrix:

$$
\begin{equation}
    \mathbf{P}=\Sigma^{-1}.
\end{equation}
$$

Given the Cholesky decomposition in Equation $21$, the Cholesky factor of the precision matrix can be written as

$$
\begin{equation}
\begin{aligned}
    \mathbf{P}&=\left(\mathbf{L}\mathbf{L}^{\top}\right)^{-1}\\
              &=\left(\mathbf{L}^{\top}\right)^{-1}\mathbf{L}^{-1}\\
              &=\left(\mathbf{L}^{-1}\right)^{\top}\mathbf{L}^{-1}\\
\end{aligned}
\end{equation}
$$

Thus, the multivariate analogue of the ‘z-score’ is

$$
\begin{equation}
    \mathbf{Z}=\mathbf{L}^{-1}(\mathbf{X}-\mathbb{E}[\mathbf{X}]).
\end{equation}
$$

Geometrically, this operation linearly transforms (rotates) the data—i.e., samples of the random vector with covariance $\Sigma$—into a new set of variables whose covariance is the identity matrix (see Figure $2$).

{% include widgets/figure.html num=2 src="/assets/images/2025-12-28/whitening.png" alt="whitening" caption="(Left) Samples of a random variable $X$ with variances $\sigma_1=1.5$ and $\sigma_2=0.5$ and correlation $\rho=0.5$. (Right) Samples rotated by the Cholesky factor of the known covariance matrix, as in Equation $26$. The effect is that the rotated samples are uncorrelated with unit variance." %}

Ignoring the mean for simplicity, we can readily verify this transformation:

$$
\begin{equation}
\begin{aligned}
\mathbb{V}[Z]&=\mathbb{E}\left[\left(\mathbf{L}^{-1}\mathbf{X}\right)\left(\mathbf{L}^{-1}\mathbf{X}\right)^\top\right]\\
        &=\mathbb{E}\left[\mathbf{L}^{-1}\mathbf{X}\mathbf{X}^\top\left(\mathbf{L}^{-1}\right)^\top\right]\\
        &=\mathbf{L}^{-1}\mathbb{E}\left[\mathbf{X}\mathbf{X}^\top\right]\left(\mathbf{L}^{-1}\right)^\top\\
        &=\mathbf{L}^{-1}\Sigma\left(\mathbf{L}^\top\right)^{-1}\\
        &=\mathbf{I}\\
\end{aligned}
\end{equation}
$$

Including the mean makes the derivation marginally more involved, but the main idea remains the same. While a full explanation of why the Cholesky decomposition works warrants a separate discussion, some important ideas are also encapsulated in **principal component analysis** (PCA).

#### Summary statistics

It is often desirable to summarise the information in a covariance matrix with a single number. To my knowledge, there are at least two common summary statistics of this kind, each reflecting different aspects.

**Total variance.** The _total variance_ of a random vector $\mathbf{X}$ is given by the trace of its covariance matrix:

$$
\begin{equation}
    \sigma_{\text{tv}}^2:=\text{tr}(\Sigma)=\sum_{i=1}^n\sigma_{i}^2.
\end{equation}
$$

This number summarises the total variance across all components of $\mathbf{X}$. The concept is important in PCA, where total variance is preserved through the transformation. Of course, in the one-dimensional case, total variance just reduces to the usual variance, that is, $\sigma_{\text{tv}}^2=\sigma^2$.

**Generalised variance.** The _generalised variance_ ([Wilks, 1932](URL)) of a random vector $\mathbf{X}$ is defined as the determinant of its covariance matrix:

$$
\begin{equation}
    \sigma_{\text{gv}}^2:=\text{det}(\Sigma)=|\Sigma|.
\end{equation}
$$

The determinant has compelling [geometric](https://www.youtube.com/watch?v=Ip3X9LOh2dk) [interpretations](), but, for our purposes, the simplest way to understand it here is as the product of the eigenvalues of $\Sigma$:

$$
\begin{equation}
    |\Sigma|=\prod_{i=1}^n\lambda_i.
\end{equation}
$$

Thus, the generalised variance reflects the overall magnitude of the linear transformation encapsulated by $\Sigma$.

Generalised variance differs sharply from total variance. For example, consider a two-dimensional covariance matrix in which

$$
\begin{equation}
    \sigma_1=\sigma_2=\rho=1
\end{equation}
$$

Here, the total variance is two, while the determinant is one. If, instead, the variables are highly correlated, say with $\rho=0.98$, then the total variance remains two, but the determinant drops to about $0.039$, indicating near-singularity. Thus, total variance simply describes the overall spread in $\Sigma$, whilst the generalised variance captures how variables in $\mathbf{X}$ co-vary; when these variables are highly (un)correlated, the generalised variance becomes low (or high).

#### Example

Let us end with two examples illustrating these ideas.

**Multivariate normal.** First, recall the probability density function (PDF) for a standard normal random variable:

$$
\begin{equation}
    p\left(x;\mu,\sigma^2\right)=\frac{1}{\sqrt{2\pi}\sigma}\exp\left\lbrace -\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2 \right\rbrace.
\end{equation}
$$

We immediately see that the squared term is simply the square from Equation $23$. Now, viewing the covariance matrix as a high-dimensional variance, the PDF for a multivariate normal random variable becomes

$$
\begin{equation}
    p\left(\mathbf{x};\mu,\Sigma \right)=\frac{1}{(2\pi)^{-n/2}|\Sigma|^{-1/2}}\exp\left\lbrace -\frac{1}{2}(\mathbf{x}-\mu)^{\top}\Sigma^{-1}(\mathbf{x}-\mu) \right\rbrace.
\end{equation}
$$

Here, the [Mahalanobis distance](https://en.wikipedia.org/wiki/Mahalanobis_distance) serves as a multivariate analogue to the whitening transformation, and the variance (the normalising constant)

$$
\begin{equation}
    \frac{1}{\sqrt{2\pi}\sigma}
\end{equation}
$$

has a multivariate form as the generalised variance.

**Correlated random variables.** Consider a scalar random variable $Z$ with unit variance, i.e. $\mathbb{V}[Z]=1$. Multiplying $Z$ by $\sigma$ yields a random variable $X$ with variance $\sigma^2$:

$$
\begin{equation}
    \mathbb{V}[\sigma Z]=\sigma^2\mathbb{V}[Z]=\sigma^2.
\end{equation}
$$

How does this extend to higher dimensions? Suppose we have two independent random variables

$$
\begin{equation}
    Z=
    \begin{bmatrix}
        Z_1 \\
        Z_2 \\
    \end{bmatrix},
\end{equation}
$$

and we wish to transform them into a random vector $\mathbf{X}$ with a specified covariance matrix $\Sigma$. The standard method is to multiply $\mathbf{Z}$ by the Cholesky factor as in Equation $26$.

This leads to a general algorithm for constructing correlated random variables: generate $\mathbf{Z}$ with independent unit variances, and multiply by the Cholesky factor of the desired covariance matrix. In the $2\times2$ case,

$$
\begin{equation}
    \text{cholesky}\left(\begin{bmatrix}1&\rho\\ \rho&1\\\end{bmatrix} \right)=\mathbf{L}\mathbf{L}^{\top}=\begin{bmatrix}1&0\\ \rho&\sqrt{1-\rho^2}\end{bmatrix} \begin{bmatrix}1&\rho\\ 0&\sqrt{1-\rho^2} \end{bmatrix}
\end{equation}
$$

Thus, the algorithm is: draw two independent (i.i.d.) random variables $Z_1$ and $Z_2$, each with unit variance, and set

$$
\begin{equation}
\begin{aligned}
    X_1&:=Z_1\\
    X_2&:=Z_1\rho+Z_2\sqrt{1-\rho^2}\\
\end{aligned}
\end{equation}
$$

This procedure conveniently generalises to any dimensionality, and non-unit variances can also be accommodated as required.

#### Conclusion

With the appropriate perspective, it is quite natural to interpret $\Sigma$ simply as ‘variance’, denoted $\mathbb{V}[\mathbf{X}]$. This view is beneficial as it renders certain properties of covariance matrices nearly self-evident, such as their positive semi-definiteness or the appearance of their inverses in whitening transformations. However, high-dimensional variance also possesses properties that are not relevant in the univariate case, such as the correlations between different variables in $\mathbf{X}$. In my opinion, the most fruitful framing is to regard univariate variance as merely a special case of covariance matrices; nevertheless, either perspective aids in developing a deeper intuition for the subject.
