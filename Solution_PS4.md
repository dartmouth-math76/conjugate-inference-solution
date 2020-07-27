---
layout: default
title: "PS3 Solution"
---

# Solution to Conjugate Inference Problem Set

## Problem 1
### Part 1:
Calculate the mean and variance of a continuous random variable $X$ with sample space $[0,1]$ and density $f(x) = 3x^2$.

#### Solution:
$$
\begin{eqnarray}
\text{Mean}[X] &=& \mathbb{E}[X]\\
&=& \int_0^1 x f(x) dx\\
&=& \int_0^1 3x^3 dx\\
& = & \left. \frac{3}{4} x^4\right|_0^1\\
& = & \frac{3}{4}
\end{eqnarray}
$$

$$
\begin{eqnarray}
\text{Var}[X] &=& \mathbb{E}[(X-\mathbb{E}[X])^2]\\
&=& \int_0^1 \left(x-\frac{3}{4}\right)^2 f(x) dx\\
&=& \int_0^1 \left(x^2-\frac{3}{2}x + \frac{9}{16}\right) f(x) dx\\
&=& \int_0^1 3x^4-\frac{9}{2}x^3 + \frac{27}{16}x^2\, dx\\
&=& \left. \left(\frac{3}{5}x^5-\frac{9}{8}x^4 + \frac{27}{48}x^3\right)\right|_0^1\\
& \approx & 0.0375 
\end{eqnarray}
$$

### Part 2:
You roll two dice. What is the expected number of sixes that will show?

#### Solution:
There are many ways of completing this task.  Here, we'll look a solution using expectations and indicator functions.   Let $X_1$ be a random variable representing the outcome of the first dice and let $X_2$ be a random variable representing the outcome of the second dice.  The sample space for both of these random variables is $\Omega = \{1,2,\ldots, 6\}$.  The number of sixes in a roll of both dice will be 
$$
g(x_1,x_2) = I[x_1=6] + I[x_2=6],
$$
where $I[\cdot]$ is an indicator function with value $1$ when the argument is true and value $0$ otherwise.   The expected number of sixes is the expected value of $g(x_1,x_2)$, which is 
$$
\mathbb{E}[g(x_1,x_2)] = \sum_{i=1}^6 \sum_{j=1}^6 g(i,j) m(i,j),
$$
where $m(x_1,x_2)$ is the joint mass function for the random variables $X_1$ and $X_2$.  Under the reasonable assumption that random variables $X_1$ and $X_2$ are independent and that $m(x_1)=\frac{1}{6}$ is uniform, we have
$$
\begin{eqnarray}
\mathbb{E}[g(x_1,x_2)] &=& \frac{1}{36}\sum_{i=1}^6 \left(\sum_{j=1}^6 \left(I[i=6] + I[j=6]\right)\right)\\
&=& \frac{1}{36}\sum_{i=1}^6 \left(6 I[i=6] + \sum_{j=1}^6I[j=6]\right)\\
&=& \frac{1}{36}\sum_{i=1}^6 \left(6 I[i=6] + 1\right)\\
&=& \frac{6+6}{36} = \frac{1}{3}
\end{eqnarray}
$$
Note that the same result could also be achieved more easily by observing that the expected number of sixes when rolling two dice is just twice the expected number of sixes when rolling one dice.  This is true by the linearity of the expectation.

### Part 3:
Consider a spring mass system where a block is attached to an immovable wall by a spring.  You apply a constant force $f$ to the mass pushing towards the left.  The mass then reaches equilibrium at a displacement $\delta x$.  At this point, the spring provides a force $K\delta x$, where the stiffness $K$ is a parameter describing the strength of the spring.  This force balances the applied force $f$ so that $f=K\delta x$.   Now consider a situation where $f$ is known but the value of $K$ follows a [Gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution).   What is the expected value of $\delta x$?

#### Solution:
First, solving for $\delta x$ we have
$$
\delta x = \frac{f}{K}
$$

##### Approach 1:
Let $g(K)$ denote the Gamma density of $K$ with paramters $\alpha$ and $\beta$ (Note the use of the $\alpha,\beta$ form of the Gamma density here, not the $k,\theta$ form we have used before).  Thus, 

$$
\begin{eqnarray}
\mathbb{E}[\delta x]&=&\mathbb{E}\left[\frac{f}{K}\right]\\
&=&\int_0^\infty \frac{f}{K} g(K) dK\\
& = & \int_0^\infty \frac{f}{K}\frac{\beta^\alpha}{\Gamma(\alpha)} K^{\alpha-1}e^{-\beta K} dK\\
&=&\frac{f\beta^\alpha}{\Gamma(\alpha)} \int_0^\infty K^{\alpha-2}e^{-\beta K} dK
\end{eqnarray}
$$

The to computing this integral is to recognize that it is similar to the [definition of the Gamma function](https://en.wikipedia.org/wiki/Gamma_function).  Setting $t=\beta K$ and performing integration by parts results in
$$
\begin{eqnarray}
\mathbb{E}[\delta x] & =& \frac{f\beta^\alpha}{\Gamma(b)} \int_0^\infty K^{\alpha-2}e^{-\beta K} dK\\
 & =& \frac{f\beta^\alpha}{\Gamma(\alpha)} \int_0^\infty  \beta^{2-\alpha} t^{\alpha-2}e^{-t} \beta^{-1}dt\\
 &=& \frac{f\beta}{\Gamma(\alpha)}\left[\left(\frac{1}{\alpha-1}t^{\alpha-1}e^{-t}\right)_0^\infty + \frac{1}{\alpha-1}\int_0^\infty t^{\alpha-1}e^{-t} dt \right]\\
 &=& \frac{f\beta}{\Gamma(\alpha)(\alpha-1)}\int_0^\infty t^{\alpha-1}e^{-t} dt\\
 &=& \frac{f\beta}{\alpha-1}
\end{eqnarray}
$$

##### Approach 2:
A much simpler alternative to the previous approach is to observe that if $K$ follows a Gamma distribution then $K^{-1}$ follows an inverse Gamma distribution, which has a known mean of $\beta/(\alpha-1)$   Thus, by the linearity of expectation
$$
\begin{eqnarray}
\mathbb{E}[\delta x] & = & \mathbb{E}\left[\frac{f}{K}\right]\\
&=& f\mathbb{E}\left[K^{-1}\right]\\
&=& \frac{f\beta}{\alpha-1}
\end{eqnarray}
$$


## Problem 2: Product of Independent RVs

Consider two independent continuous random variables $X$ and $Y$ with sample spaces $\Omega_x$, $\Omega_y$, and densities $f(x)$, $f(y)$.   Show that the expectation of the product $XY$ is the product of the expectations $\mathbb{E}[X]\mathbb{E}[Y]$.

#### Solution:
By the definition of independence, we know the joint density $f(x,y)=f(x)f(y)$.    Thus
$$
\begin{eqnarray}
\mathbb{E}[XY] & = & \int_{\Omega_y}\int_{\Omega_x} xy f(x,y) dx dy\\
 & = & \int_{\Omega_y}\int_{\Omega_x} xy f(x) f(y) dx dy\\
 & = & \int_{\Omega_y} y f(y) \int_{\Omega_x} x f(x) dx dy\\
 & = & \left( \int_{\Omega_y} y f(y) dy\right)\left( \int_{\Omega_x} x f(x) dx\right)\\
 &=& \mathbb{E}[Y]\mathbb{E}[X]
\end{eqnarray}
$$


## Problem 3: NYC Bicycling Trends

On average, are more bikes likely to cross the Brooklyn Bridge in May or September?  Answering that question is the the goal of this problem.   New York City shares a lot of their data on the [OpenData](https://opendata.cityofnewyork.us/) website.  Here, we will analyze data from 2017 that counts the number of bikers crossing NYC bridges each day from April through October.  

The data collected by NYC counts the number of bikers crossing bridges each day.  To model this in our likelihood function, we will employ the Poisson distribution, which is commonly used to describes counts per day or counts per area.  The Poisson distribution describes the probability mass function for a discrete random variable that has a sample space which includes all positive integers.   Let $B$ be a Poisson random variable representing the number of bicycles crossing the Brooklyn Bridge.   The Poisson mass function $m(k)$ takes the form

$$
m(b|\lambda) = \mathbb{P}[B=b] = \frac{\lambda^b e^{-\lambda}}{b!}
$$

The Poisson distribution has a single parameter $\lambda$ that represents the average rate (e.g., bicyclists per day) of the random variable $B$.  Interestingly, both the mean and variance of $B$ are equal to $\lambda$.  Our goal here is to define a posterior density $f(\lambda | b_{May})$ using observations of the bicyclist count in May,  i.e., $b_{May}$, and another posterior density $f(\lambda | b_{Sep})$ using observations of the bicyclist count in September.  Comparing the two posteriors will allow us to compute the probability that more bikers will cross Brooklyn Bridge in May than September.

We will assume that the prior density $f(\lambda)$ is a Gamma density.  The Gamma distribution is a flexible family of distributions useful for describing positive real-valued variables (like $\lambda$).  The probability density function of a Gamma random variable is given by

$$
f(\lambda) = \text{Gamma}(\lambda; k, \theta) = \frac{1}{\Gamma(k)\theta^k} \lambda^{k-1} e^{-\frac{\lambda}{\theta}},
$$

where $k$ and $\theta$ are parameters describing the shape of the density.  

### Part 1:
To show that the Gamma distribution is conjugate for the Poisson likelihood, we can write out the form of the posterior density and look for connections with the Gamma density.

$$
\begin{eqnarray}
f(\lambda | b) &=& \frac{m(b|\lambda) f(\lambda)}{m(b)}\\
&=& C_1 m(b|\lambda) f(\lambda)\\
& = & C_1 \frac{\lambda^b e^{-\lambda}}{b!} \frac{1}{\Gamma(k)\theta^k} \lambda^{k-1} e^{-\frac{\lambda}{\theta}}\\
& = & C_2 \lambda^be^{-\lambda}\lambda^{k-1} e^{-\frac{\lambda}{\theta}}\\
& = & C_2 \lambda^{b+k-1}e^{-\lambda-\frac{\lambda}{\theta}}\\
& = & C_2 \lambda^{k^\ast-1}e^{-\frac{\lambda}{\theta^\ast}},
\end{eqnarray}
$$

where $k^\ast = b+k$ and 

$$
\begin{eqnarray}
-\frac{\lambda}{\theta^\ast} &= & -\lambda-\frac{\lambda}{\theta} \\
\Rightarrow \theta^\ast &=& \frac{\lambda}{\lambda + \frac{\lambda}{\theta}}\\
&=& \frac{\lambda \theta}{\lambda\theta + \lambda}\\
&= & \frac{\theta}{\theta+1}.
\end{eqnarray}
$$

### Part 2:
Using the result from part 1, derive the parameters $k^\ast,\theta^\ast$ of the Gamma *posterior* given $N$ independent observations $b_1,b_2,\ldots, b_N$, where $b_i$ denotes the number of bicyclists that crossed Brooklyn Bridge on day $i$.  Recall that we are assuming the likelihood $f(b_i | \lambda)$ is a Poisson mass function. 

*Hint:* write $k^\ast$ and $\theta^\ast$ in terms of $b_1,\ldots,b_N$, $k$, and $\theta$.

#### Solution:
Let $k_i$ and $\theta_i$ denote the parameters of the posterior Gamma density after $i$ observations.   From the previous part of this problem, we have  the following recursive relationships:
$$
k_1 = b_1+k_{0}
$$
and 
$$
\theta_1 = \frac{\theta_{0}}{\theta_{0} + 1},
$$
where, in general, $b_i$ is the value of observation $i$ (i.e., the number  of bikers observed on day $i$).

For $k_N$, we would therefore have 
$$
k_N=k_0 +\sum_{i=1}^Nb_i.
$$

The relationship for $\theta$ is less obvious however. Looking at another step of the $\theta$ update, we have 
$$
\begin{eqnarray}
\theta_{2} &=& \frac{\theta_{1}}{\theta_{1} + 1},\\
&=& \frac{\frac{\theta_{0}}{\theta_{0} + 1}}{\frac{\theta_{0}}{\theta_{0} + 1} + 1}\\
&=&  \frac{\theta_{0}}{2\theta_{0} + 1}
\end{eqnarray}
$$

and again
$$
\begin{eqnarray}
\theta_{3} &=& \frac{\theta_{2}}{\theta_{2} + 1},\\
&=& \frac{\frac{\theta_{0}}{2\theta_{0} + 1}}{\frac{\theta_{0}}{2\theta_{0} + 1} + 1}\\
&=&  \frac{\theta_{0}}{3\theta_{0} + 1}.
\end{eqnarray}
$$

In general, we therefore have
$$
\begin{eqnarray}
\theta_{N} &=& \frac{\theta_{0}}{N\theta_{0} + 1}
\end{eqnarray}
$$



### Part 3:
Read through the code below to and use the result of part 2 to update the `PosteriorParameters` function.  This is labeled as `TODO` in the code.

#### Solution:
See code in [this repository](https://github.com/dartmouth-math76/conjugate-inference-solution).

### Part 4:
Using your code from Part 3, run the code below that computes the posterior parameters $k^\ast_{May},\theta^\ast_{May}$ and $k^\ast_{Sep},\theta^\ast_{Sep}$ using data from May and September, respectively.  Plots of the posterior densities will also be generated. 

Once you've found the posterior parameters, compute the probability that a random variable $\lambda_{May}$ with density $\text{Gamma}(\lambda; k^\ast_{May}, \theta^\ast_{May})$ is greater than a random variable $\lambda_{Sep}$ with density $\text{Gamma}(\lambda; k^\ast_{Sep}, \theta^\ast_{Sep})$.  

*Hint:*  For two independent Gamma random variables $X$ and $Y$ with the same scaling parameter $\theta$ but different parameters $k_x$ and $k_y$, the random variable $Z=X/(X+Y)$ is distributed according to a [Beta distribution](https://en.wikipedia.org/wiki/Beta_distribution) with parameters $k_x$ and $k_y$.  More precisely,

$$
\frac{X}{X+Y} = Z \sim \text{Beta}(k_x,k_y)
$$

#### Solution:
From Part 2 of this problem, we see that the posterior value of $\theta$ depends only on the *number* of observations not on the actual  value of the observations.   Because May and September have the same number of days, we will have the same number of observations for both months, and the posterior values $\theta_{May}$ and $\theta_{Sep}$ will be the  same.   It is therefore possible to use the identity above.

Consider the ratio 
$$
Z =\frac{\lambda_{Sep}}{\lambda_{Sep}+\lambda_{May}}
$$
If $Z\leq \frac{1}{2}$, then $\lambda_{May}\geq\lambda_{Sep}$.  From the identity described in the project description, we know that the ratio $Z$ is distributed according to a Beta distribution with parameters $a=k_{Sep}$ and $b=k_{May}$.   Let $F(z; a,b)$ denote the CDF of a Beta random variable with parameters $a$ and $b$.  Then we have 
$$
\begin{eqnarray}
\mathbb{P}[\lambda_{May}\geq\lambda_{Sep}] & = & \mathbb{P}[Z\leq\frac{1}{2}]\\
& = & F\left(\frac{1}{2};\lambda_{Sep}, \lambda_{May}\right)
\end{eqnarray}
$$

This expression is evaluated in the code to compute the probability that the expected bicycle crossings in May $\lambda_{May}$ is greater than the expected number of bicycle crossings in September $\lambda_{Sep}$.   Using the code, we see that the probability is very small:
$$
\mathbb{P}[\lambda_{May}\geq\lambda_{Sep}] \approx 0.0046
$$

This makes sense when looking at the posterior probability density functions; the bulk of the September posterior density function is shifted to the right of the May psoterior density function.

### Part 5:
Do you think your answer from Part 4 answers the question: "What is the probability that more bikers will cross the Brooklyn Bridge in May then September?"  Why or why not?

#### Solution:
The parameter $\lambda$ in the Poisson distribution represents the. average number of bicyclists expected per day.  It is the average rate.    Our result from part 4 therefore answers  the questions "What is the probability that the average number  of bikers crossing Brooklyn Bridge per day in May will be larger than the average number of bikers crossing Brooklyn Bridge per day in September?"    This result about the rate is slightly different than asking for the probability that more bikers will cross the bridge.  The Poisson likelihood function describes the distribution over bicycle crossings on a single day for a fixed parameter $\lambda$.   Of course, we have a distribution over $\lambda$ as well.   To characterize the distribution over bikers per day, we would therefore need to compute the following marginal densities  e.g., 
$$
m(b_i | \text{May}) = \int_0^\infty m(b_i | \lambda_{May}) f(\lambda_{May})
$$
$$
m(b_i | \text{Sep}) = \int_0^\infty m(b_i | \lambda_{Sep}) f(\lambda_{Sep})
$$


