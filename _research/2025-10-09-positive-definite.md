---
title: Understanding Positive Definite Matrices
date: 2025-10-09
layout: blog
cover: /assets/images/2025-10-09/definiteness.png
tags:
  - Linear Algebra
intro: A discussion on geometric interpretation of positive definite matrices and how this relates to various properities of them, such as positive eigenvalues, positives determinants, and decomposability. Also a brief discussion on quadratic programming.
---

A real-valued matrix $\textbf{A}$ is _positive definite_ if, for every real-valued vector $\textbf{x}$,

$$
\begin{equation}
    \textbf{a}^\top\textbf{Ax}>0,\quad\textbf{x}\neq 0
\end{equation}
$$

The matrix is _positive semi-definite if_

$$
\begin{equation}
    \textbf{x}^\top\textbf{Ax}\geq 0.
\end{equation}
$$

If the inequalities are reversed, then $\textbf{A}$ is _negative definite_ or _negative semi-definite_ respectively. If no inequality holds, the matrix is _indefinite_. These definitions can be generalised to complex-valued matrices, but we will work with real numbers.

The notion of positive definiteness comes up a lot in statistics, where covariance matrices are always positive semi-definite. However, it is difficult to find a real institution for this definition. What do these definitions represent geometrically? What constraints imply that covariance matrices are always positive semi-definite? Why are these matrices featured heavily in quadratic programming? Why is the space of all positive semi-definite matrices a cone? The goal of this post is to fix this deficit in understanding.

This post will rely heavily on a geometrical understanding and matrices and dot products.

#### Multivariate positive numbers

To start, let us put aside the definition of positive semi-definite matrices and instead focus on an intuitive idea of what they are.

Positive definite matrices can be viewed as multivariate analogs to strictly positive real numbers, while positive semi-definite matrices can be viewed as multivariate analogs to nonnegative real numbers. Consider a scalar $a$. If $a<0$, then the sign of $ab$ will depend on the sign of $b$. However, if $a>0$, then $ab$ will have the same sign as $b$. We can think about $a$ and $b$ as $1$-vectors and about the product $ab$ as a dot product,

$$
\begin{equation}
    \textbf{a}^\top\textbf{b}=\large[a\large]^\top\large[b\large]=ab.
\end{equation}
$$

When $\textbf{a}$ is positive, then $\textbf{a}^\top\textbf{b}$ will be on the same side of the origin $0$ as $\textbf{b}$. So $\textbf{a}$ does not "flip" $\textbf{b}$ about the origin. Rather, it stretches $\textbf{b}$ in the same direction (Figure $1$).

{% include widgets/figure.html num=1 src="/assets/images/2025-10-09/number_line.png" alt="number line" caption="When $a>0$, the product $ab$ can be interpreted as $b$ being stretched by $a$ in the same direction that $b$ already points, where $a,b$ and $ab$ are viewed as vectors with tails at the origin $0$. On the $y$-axis, vector are offset from $0$ for visualisation purposes." %}

Analogously, a positive definite matrix behaves like a positive number in the sense that it never flips a vector about the origin $0$. The simplest example of a positive definite matrix is a diagonal matrix that scales a vector in the direction that it already points, and the simplest example of a matrix that is not positive definite is one that reverses the vector (Figure $2$).

{% include widgets/figure.html num=2 src="/assets/images/2025-10-09/obvious_pd.png" alt="pd" caption="(Left) A diagonal matrix $\textbf{A}$ with positive diagonal values does not flip $\textbf{x}$ about the origin and is thus positive definite. (Right) A diagonal matrix $\textbf{A}$ with negative values flips $\textbf{x}$ about the origin and is thus not positive definite." %}

However, it is less clear how to think about scenarios in which a matrix $\textbf{A}$ projects in a way that does not simply scale the vector by some constant. For example, does the matrix in the left subplot of Figure $3$ flip the vector $\textbf{x}$ about the origin or not? What about the matrix in the middle subplot? How can we make these intuitions precise? The answer is in the dot product. The _dot product_ between two $N$-vectors, defined as

$$
\begin{equation}
    \mathbf{x}^\top \mathbf{y}=\sum_{i=1}^{N}x_ny_n,
\end{equation}
$$

can be viewed as vector-matrix multiplication, which in turn can be viewed as a linear production: an $N$-vector $\mathbf{y}$ is projected onto a real number line defined by the columns of a $1\times N$ matrix $m$. Another definition of the dot product is that it is the product of the Euclidean magnitudes of the two vectors times the cosine of the angle between them:

$$
\begin{equation}
    \mathbf{x}^\top \mathbf{y}=\Vert\mathbf{x}\Vert \mathbf{y}\Vert\cos\theta.
\end{equation}
$$

What these definitions imply is that two nonzero vectors $\mathbf{x}$ and $\mathbf{y}$ are _orthogonal_ to each other if their dot product is zero,

$$
\begin{equation}
    \mathbf{x}^\top \mathbf{y}=0.
\end{equation}
$$

Why? First, if their dot product is zero and if neither vector is the zero vector, then the cosine angle between the two vectors must be $90\degree$ (Figure $3$, middle).

{% include widgets/figure.html num=3 src="/assets/images/2025-10-09/ambiguous_pd.png" alt="ambiguous" caption="A vector $\mathbf{x}$ is mapped by a matrix $\mathbf{A}$ to a new vector $\mathbf{Ax}$. (Left) The angle between $\mathbf{x}$ and $\mathbf{Ax}$ is acute. (Middle) The angle between $\mathbf{x}$ and $\mathbf{Ax}$ is right. The angle between $\mathbf{x}$ and $\mathbf{Ax}$ is obtuse." %}

If the dot product between two nonzero vectors is positive.

$$
\begin{equation}
    \mathbf{x}^\top \mathbf{y}>0,
\end{equation}
$$

then both vectors point in a "similar direction" in the sense the cosine angle between then is less than $90\degree$; the angle is acute (Figure 3, left). And if the dot product is negative,

$$
\begin{equation}
    \mathbf{x}^\top \mathbf{y}<0,
\end{equation}
$$

then both vectors point in a "dissimilar direction" in the sense that the cosine angle between them is greater then $90\degree$; the angle is obtuse (Figure 3, right).

As we can see, one geometric interpretation of positive definiteness is that $\mathbf{A}$ does not "flip" $\mathbf{x}$ about the origin in the sense that the cosine angle between $\mathbf{x}$ and $\mathbf{Ax}$ is between $-90\degree$ and $90\degree$. In other words, $\mathbf{x}$ and $\mathbf{Ax}$ are in the same half of an $N$-dimensional hypersphere (Figure 4).

{% include widgets/figure.html num=4 src="/assets/images/2025-10-09/cosine_angle.png" alt="cosine_angle" caption="(Left) The cosine function is positive between $-\pi/2$ and $-\pi/2$ or between $-90\degree$ and $90\degree$. (Right) A positive semi-definite matrix $\mathbf{A}$ maps a vector $\mathbf{x}$ to a new location $\mathbf{Ax}$ such that the cosine angle between $\mathbf{x}$ and $\mathbf{Ax}$ is between $-90\degree$ and $90\degree$." %}

Figure $4$ captures the geometric interpretation of Equation $1$. It is a easier to see this when we make the location of the dot product operator more explicit, i.e.

$$
\begin{equation}
    \mathbf{x}\cdot \mathbf{Ax}>0.
\end{equation}
$$

(Here, $\mathbf{a}\cdot \mathbf{b}$ is another common notation for the dot product.) This is similar to defining a real number as a number $a$ such that

$$
\begin{equation}
    b(ab)>0,
\end{equation}
$$

for any number $b\neq 0$. To summarise everything so far, an $N\times N$ matrix is positive definite if it maps any vector into the same half of an $N$-dimensional hypersphere. Geometrically, this is analogous to a positive number (Figure $1$).

#### Examples

How that we have a definition, as well as an intuition, for positive semi-definite matrices, let's look at some examples. First, consider the matrix

$$
\begin{equation}
\mathbf{A}=
    \begin{bmatrix}
    1 & 0 \\
    0 & 3
    \end{bmatrix}
\end{equation}
$$

This is a positive semi-definite sine
