#' Learning the causal structure of a nonlinear ZiG-DAG
#' 
#' \code{nonlinear.zidag} learns the causal structure of a nonlinear ZiG-DAG using a hill-climbing greedy search algorithm. Specifically, the conditional distribution of each node in the model is assumed to be a zero-inflated hyper-Poisson distribution. 
#'
#' @param dat a data matrix.
#' @param start a square adjacency matrix, the directed acyclic graph to be used to initialize the algorithm. If none is specified, the empty graph is used. 
#' @param nbasis an integer, the number of cubic B-spline basis functions to be used for the nonlinear ZiG-DAG. 
#' @param maxiter an integer, the maximum number of iterations.
#' @param tol a numeric value, the tolerance for the convergence of the network score (BIC).
#' @param optim.control a list of control parameters passed to \code{optim}.
#' @param verbose a boolean value. If \code{TRUE}, progress of the algorithm is printed; otherwise the function is completely silent. 
#'
#' @return An object of class \code{nonlinear.zidag}, a list containing the following components:\itemize{
#' \item\code{est}: a list of model parameter estimates for the nonlinear ZiG-DAG, of which the component "\code{E}" gives the adjacency matrix of the estimated DAG for the nonlinear ZiG-DAG.\cr
#' \item\code{bic}: the Bayesian Information Criterion for the estimated nonlinear ZiG-DAG model.
#' \item\code{iter}: the number of iterations of the hill-climbing greedy search algorithm used. 
#' @export
#'
#' @examples
nonlinear.zidag = function(dat, start, nbasis = 4, maxiter = 500, tol = .Machine$double.eps^0.25, optim.control = list(), verbose = FALSE)
{
   n = nrow(dat) 
   p = ncol(dat)
   
   control = list(fnscale = -1, maxit = 10000, reltol = 1.0e-8)
   control[names(optim.control)] = optim.control
   
   if (is.null(start))
   {
      start = matrix(0, p, p)
   }
   
   # generate the B-spline basis matrix
   B = array(0, dim = c(n, p, nbasis))
   for (j in 1 : p)
   {
      B[ , j, ] = splines::bs(dat[ , j], df = nbasis)
   }
   
   # fit the non-linear ZIHPBN given the starting DAG
   bic_curr = rep(NA, p)
   est_curr = list()
   est_curr$E      = starting.dag
   est_curr$phi    = array(0, dim = c(p, p, nbasis))
   est_curr$psi    = array(0, dim = c(p, p, nbasis))
   est_curr$delta  = rep(0, p)
   est_curr$gamma  = rep(0, p)
   est_curr$lambda = rep(1, p)
   
   for (j in 1 : p)
   {
      pa_j    = (est_curr$E[j, ] == 1) 
      p_j     = sum(pa_j)
      start_j = c(t(est_curr$phi[j, pa_j, ]), est_curr$delta[j], t(est_curr$psi[j, pa_j, ]), est_curr$gamma[j], log(est_curr$lambda[j]))
      if (p_j > 0) 
      {
         X_pa_j = t(apply(B[ , pa_j, ], 1, c))
      } else 
      {
         X_pa_j = matrix(NA, n, 0)
      }
      
      out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = X_pa_j, lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                      gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = X_pa_j), method = "BFGS", control = control)
      bic_curr[j] = -2 * out_j$value + (2 * p_j * nbasis + 3) * log(n)
      
      if (p_j > 0)
      {
         est_curr$phi[j, pa_j, ] = matrix(out_j$par[1 : (p_j * nbasis)], p_j, nbasis, byrow = TRUE)
         est_curr$psi[j, pa_j, ] = matrix(out_j$par[(p_j * nbasis + 2) : (2 * p_j * nbasis + 1)], p_j, nbasis, byrow = TRUE) 
      }
      est_curr$delta[j]  = out_j$par[p_j * nbasis + 1]
      est_curr$gamma[j]  = out_j$par[2 * p_j * nbasis + 2]
      est_curr$lambda[j] = exp(out_j$par[2 * p_j * nbasis + 3])
   }
   
   # start the hill climbing
   bic_iter = bic_curr
   est_iter = est_curr
   
   for (iter in 1 : maxiter)
   {
      IMPROV  = FALSE
      
      # search all DAGs reachable from the current DAG
      # addition or deletion of an edge 
      for (j in 1 : p)
      {
         for (k in 1 : p)
         {
            if (j == k) next
            
            if (est_curr$E[j, k] == 0)
            {
               E_cand = est_curr$E
               E_cand[j , k] = 1
               G_cand = igraph::graph_from_adjacency_matrix(E_cand)
               if (!igraph::is.dag(G_cand)) next
               
               pa_j    = (E_cand[j, ] == 1) 
               p_j     = sum(pa_j)
               start_j = c(t(est_curr$phi[j, pa_j, ]), est_curr$delta[j], t(est_curr$psi[j, pa_j, ]), est_curr$gamma[j], log(est_curr$lambda[j]))
               if (p_j > 0) 
               {
                  X_pa_j = t(apply(B[ , pa_j, ], 1, c))
               } else 
               {
                  X_pa_j = matrix(NA, n, 0)
               }
               
               out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = X_pa_j, lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                               gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = X_pa_j), method = "BFGS", control = control)
               bic_cand = bic_curr
               bic_cand[j] = -2 * out_j$value + (2 * p_j * nbasis + 3) * log(n)
               
               if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
               {
                  IMPROV  = TRUE
                  bic_iter = bic_cand
                  est_iter = est_curr
                  est_iter$E = E_cand
                  est_iter$phi[j, pa_j, ] = matrix(out_j$par[1 : (p_j * nbasis)], p_j, nbasis, byrow = TRUE)
                  est_iter$psi[j, pa_j, ] = matrix(out_j$par[(p_j * nbasis + 2) : (2 * p_j * nbasis + 1)], p_j, nbasis, byrow = TRUE) 
                  est_iter$delta[j]  = out_j$par[p_j * nbasis + 1]
                  est_iter$gamma[j]  = out_j$par[2 * p_j * nbasis + 2]
                  est_iter$lambda[j] = exp(out_j$par[2 * p_j * nbasis + 3])
               }
            } else
            {
               E_cand = est_curr$E
               E_cand[j, k] = 0
               
               pa_j    = (E_cand[j, ] == 1) 
               p_j     = sum(pa_j)
               start_j = c(t(est_curr$phi[j, pa_j, ]), est_curr$delta[j], t(est_curr$psi[j, pa_j, ]), est_curr$gamma[j], log(est_curr$lambda[j]))
               if (p_j > 0) 
               {
                  X_pa_j = t(apply(B[ , pa_j, ], 1, c))
               } else 
               {
                  X_pa_j = matrix(NA, n, 0)
               }
               
               out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = X_pa_j, lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                               gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = X_pa_j), method = "BFGS", control = control)
               bic_cand = bic_curr
               bic_cand[j] = -2 * out_j$value + (2 * p_j * nbasis + 3) * log(n)
               
               if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
               {
                  IMPROV  = TRUE
                  bic_iter = bic_cand
                  est_iter = est_curr
                  est_iter$E = E_cand
                  est_iter$phi[j, k, ] = est_iter$psi[j, k, ] = 0
                  if (p_j > 0)
                  {
                     est_iter$phi[j, pa_j, ] = matrix(out_j$par[1 : (p_j * nbasis)], p_j, nbasis, byrow = TRUE)
                     est_iter$psi[j, pa_j, ] = matrix(out_j$par[(p_j * nbasis + 2) : (2 * p_j * nbasis + 1)], p_j, nbasis, byrow = TRUE) 
                  }
                  est_iter$delta[j]  = out_j$par[p_j * nbasis + 1]
                  est_iter$gamma[j]  = out_j$par[2 * p_j * nbasis + 2]
                  est_iter$lambda[j] = exp(out_j$par[2 * p_j * nbasis + 3])
               }
            }
         }
      }
      
      # reversal of an edge
      id_edge = which(est_curr$E == 1, arr.ind = TRUE)
      n_rev    = nrow(id_edge)
      
      if (n_rev > 0)
      {
         for (l in 1 : n_rev)
         {
            j = id_edge[l, 1]
            k = id_edge[l, 2]
            
            E_cand = est_curr$E
            E_cand[j , k] = 0
            E_cand[k , j] = 1
            G_cand = igraph::graph_from_adjacency_matrix(E_cand)
            if (!igraph::is.dag(G_cand)) next
            
            pa_j    = (E_cand[j, ] == 1) 
            p_j     = sum(pa_j)
            start_j = c(t(est_curr$phi[j, pa_j, ]), est_curr$delta[j], t(est_curr$psi[j, pa_j, ]), est_curr$gamma[j], log(est_curr$lambda[j]))
            if (p_j > 0) 
            {
               X_pa_j = t(apply(B[ , pa_j, ], 1, c))
            } else 
            {
               X_pa_j = matrix(NA, n, 0)
            }
            
            pa_k    = (E_cand[k, ] == 1) 
            p_k     = sum(pa_k)
            start_k = c(t(est_curr$phi[k, pa_k, ]), est_curr$delta[k], t(est_curr$psi[k, pa_k, ]), est_curr$gamma[k], log(est_curr$lambda[k]))
            if (p_k > 0) 
            {
               X_pa_k = t(apply(B[ , pa_k, ], 1, c))
            } else 
            {
               X_pa_k = matrix(NA, n, 0)
            }
            
            out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = X_pa_j, lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                            gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = X_pa_j), method = "BFGS", control = control)
            out_k   = optim(par = start_k, fn = function(z) dZIHP_cpp(z, y = dat[, k], X = X_pa_k, lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                            gr = function(z) gradZIHP_cpp(z, y = dat[, k], X = X_pa_k), method = "BFGS", control = control)
            bic_cand = bic_curr
            bic_cand[j] = -2 * out_j$value + (2 * p_j * nbasis + 3) * log(n)
            bic_cand[k] = -2 * out_k$value + (2 * p_k * nbasis + 3) * log(n)
            
            if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
            {
               IMPROV  = TRUE
               bic_iter = bic_cand
               est_iter = est_curr
               est_iter$E = E_cand
               est_iter$phi[j, k, ] = est_iter$psi[j, k, ] = 0
               if (p_j > 0)
               {
                  est_iter$phi[j, pa_j, ] = matrix(out_j$par[1 : (p_j * nbasis)], p_j, nbasis, byrow = TRUE)
                  est_iter$psi[j, pa_j, ] = matrix(out_j$par[(p_j * nbasis + 2) : (2 * p_j * nbasis + 1)], p_j, nbasis, byrow = TRUE)
               }
               est_iter$delta[j]  = out_j$par[p_j * nbasis + 1]
               est_iter$gamma[j]  = out_j$par[2 *p_j * nbasis + 2]
               est_iter$lambda[j] = exp(out_j$par[2 * p_j * nbasis + 3])
               est_iter$phi[k, pa_k, ] = matrix(out_k$par[1 : (p_k * nbasis)], p_k, nbasis, byrow = TRUE)
               est_iter$psi[k, pa_k, ] = matrix(out_k$par[(p_k * nbasis + 2) : (2 * p_k * nbasis + 1)], p_k, nbasis, byrow = TRUE)
               est_iter$delta[k]  = out_k$par[p_k * nbasis + 1]
               est_iter$gamma[k]  = out_k$par[2 * p_k * nbasis + 2]
               est_iter$lambda[k] = exp(out_k$par[2 * p_k * nbasis + 3])
            }
         }
      }
      
      if (verbose) 
         cat("iter =", iter, "; BIC =", round(sum(bic_iter), 4), "\n")
      
      if (!IMPROV) break
      
      bic_curr = bic_iter
      est_curr = est_iter
   }
   
   out = list()
   out$est  = est_curr
   out$bic  = sum(bic_curr)
   out$iter = iter
   clas(out) = "nonlinear.zigdag"
   return(out)
}
