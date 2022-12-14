% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonlinear.zigdag.R
\name{nonlinear.zigdag}
\alias{nonlinear.zigdag}
\title{Learning the causal structure of a nonlinear ZiG-DAG}
\usage{
nonlinear.zigdag(
  dat,
  start = NULL,
  nbasis = 4,
  maxiter = 500,
  tol = .Machine$double.eps^0.25,
  optim.control = list(),
  verbose = FALSE
)
}
\arguments{
\item{dat}{a data matrix.}

\item{start}{a square adjacency matrix, the directed acyclic graph to be used to initialize the algorithm. If none is specified, the empty graph is used.}

\item{nbasis}{an integer, the number of cubic B-spline basis functions to be used for the nonlinear ZiG-DAG.}

\item{maxiter}{an integer, the maximum number of iterations.}

\item{tol}{a numeric value, the tolerance for the convergence of the network score (BIC).}

\item{optim.control}{a list of control parameters passed to \code{optim}.}

\item{verbose}{a boolean value. If \code{TRUE}, progress of the algorithm is printed; otherwise the function is completely silent.}
}
\value{
An object of class \code{nonlinear.zidag}, a list containing the following components:\itemize{
\item\code{est}: a list of model parameter estimates for the nonlinear ZiG-DAG, of which the component "\code{E}" gives the adjacency matrix of the estimated DAG for the nonlinear ZiG-DAG.\cr
\item\code{bic}: the Bayesian Information Criterion for the estimated nonlinear ZiG-DAG model.
\item\code{iter}: the number of iterations of the hill-climbing greedy search algorithm used. 
}
}
\description{
\code{nonlinear.zigdag} learns the causal structure of a nonlinear ZiG-DAG using a hill-climbing greedy search algorithm. Specifically, the conditional distribution of each node in the model is assumed to be a zero-inflated hyper-Poisson distribution.
}
\examples{
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
}
