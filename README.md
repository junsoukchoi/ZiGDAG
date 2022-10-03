ZiGDAG: R package for model-based causal discovery for zero-inflated
count data
================

-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#example-of-linear-zig-dags"
    id="toc-example-of-linear-zig-dags">Example of linear ZiG-DAGs</a>
-   <a href="#example-of-nonlinear-zig-dags"
    id="toc-example-of-nonlinear-zig-dags">Example of nonlinear ZiG-DAGs</a>

The R package `ZiGDAG` implements zero-inflated generalized
hypergeometric directed acyclic graphs (ZiG-DAGs) for inference of
causal structure from observational zero-inflated count data. For the
structure learning of ZiG-DAGs, score-based greedy search algorithms are
implemented.

## Installation

To install the latest version from Github, use

``` r
library(devtools)
devtools::install_github("junsoukchoi/ZiGDAG")
```

## Example of linear ZiG-DAGs

``` r
library(ZiGDAG)
library(igraph)
library(DGLMExtPois)
# set a random seed
set.seed(1)

# generate synthetic data from a linear ZiG-DAG, where zero-inflated hyper-Poisson distributions are assumed for each node
# generate the DAG with 5 nodes: 1 -> 2 -> 3 -> 4 -> 5
p = 5   # number of variables
E_true = matrix(0, p, p)
for (j in 2 : p) E_true[j, j - 1] = 1

# generate model parameters for the linear ZiG-DAG
alpha_true  = matrix(0, p, p)
alpha_true[E_true == 1] = runif(sum(E_true), 0.5, 2)
beta_true   = matrix(0, p, p)
beta_true[E_true == 1]  = runif(sum(E_true), -2, -0.5)
delta_true  = runif(p, -1.5, -1)
gamma_true  = runif(p, 1, 1.5)
lambda_true = exp(runif(p, -2, 2))

# generate synthetic data from the specified linear ZiG-DAG with sample size n = 500 
n   = 500  # sample size
dat = matrix(0, n, p)
order_nodes = as_ids(topo_sort(graph_from_adjacency_matrix(t(E_true))))
for (j in order_nodes)
{
   pi_true = as.vector(exp(dat %*% alpha_true[j, ] + delta_true[j]))
   pi_true = pi_true / (1 + pi_true)
   pi_true[is.nan(pi_true)] = 1
   mu_true = as.vector(exp(dat %*% beta_true[j, ] + gamma_true[j]))
   dat[ , j] = (1 - rbinom(n, 1, pi_true)) * rhP(n, lambda_true[j], mu_true)
}

# learn the causal structure of linear ZiG-DAG from synthetic data 
# apply the hill-climbing algorithm for linear ZiG-DAG to synthetic data
fit = linear.zigdag(dat, ghpd = "hyper.poisson", method = "hc")
# adjacency matrix of the estimated causal DAG
E_est = fit$est$E
```

## Example of nonlinear ZiG-DAGs

``` r
library(ZiGDAG)
library(igraph)
library(DGLMExtPois)
# set a random seed
set.seed(1)

# generate synthetic data from a nonlinear ZiG-DAG, where zero-inflated hyper-Poisson distributions are assumed for each node
# generate the DAG with 5 nodes: 1 -> 2 -> 3 -> 4 -> 5
p = 5   # number of variables
E_true = matrix(0, p, p)
for (j in 2 : p) E_true[j, j - 1] = 1

# generate nonlinear functions and parameters for the nonlinear ZiG-DAG 
f_nonlinear = function(z, type)
{
   if (type == 1) return(0.5 * z * (z - 3))
   if (type == 2) return(sin(z))
}

g_nonlinear = function(z, type)
{
   if (type == 1) return(-0.5 * (z - 1.5)^2)
   if (type == 2) return(cos(z))
}

f_true = matrix(0, p, p)
f_true[E_true == 1] = sample(1 : 2, size = sum(E_true), replace = TRUE)
g_true = matrix(0, p, p)
g_true[E_true == 1] = sample(1 : 2, size = sum(E_true), replace = TRUE)
delta_true  = runif(p, -1.5, -1)
gamma_true  = runif(p, 1, 1.5)
lambda_true = exp(runif(p, -2, 2))

# generate synthetic data from the specified nonlinear ZiG-DAG with sample size n = 500 
n   = 500  # sample size
dat = matrix(0, n, p)
order_nodes = as_ids(topo_sort(graph_from_adjacency_matrix(t(E_true))))
for (j in order_nodes)
{
   p_j  = sum(E_true[j, ] == 1)
   pi_true = rep(delta_true[j], n)
   mu_true = rep(gamma_true[j], n)
   if (p_j > 0)
   {
      pa_j = which(E_true[j, ] == 1)
      for (a in 1 : p_j)
      {
         k = pa_j[a] 
         pi_true = pi_true + f_nonlinear(dat[ , k], type = f_true[j, k])
         mu_true = mu_true + g_nonlinear(dat[ , k], type = g_true[j, k])
      }
   }
   
   pi_true = exp(pi_true) / (1 + exp(pi_true))
   pi_true[is.nan(pi_true)] = 1
   mu_true = exp(mu_true)
   dat[ , j] = (1 - rbinom(n, 1, pi_true)) * rhP(n, lambda_true[j], mu_true)
}

# learn the causal structure of nonlinear ZiG-DAG from synthetic data 
# apply the hill-climbing algorithm for nonlinear ZiG-DAG to synthetic data
# 4 cubic B-spline functions are used to approximate nonlinear functions in nonlinear ZiG-DAG models.
fit = nonlinear.zigdag(dat, nbasis = 4)   
# adjacency matrix of the estimated causal DAG
E_est = fit$est$E
```
