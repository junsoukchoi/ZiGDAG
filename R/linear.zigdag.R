#' Learning the causal structure of a linear ZiG-DAG
#' 
#' \code{linear.zidag} learns the causal structure of a linear ZiG-DAG using a score-based greedy search algorithm.
#'
#' @param dat a data matrix.
#' @param start a square adjacency matrix, the directed acyclic graph to be used to initialize the algorithm. If none is specified, the empty graph is used. 
#' @param ghpd a character string, the generalized hypergeometric probability distribution to be used for the linear ZiG-DAG. Possible values are "\code{hyper.poisson}" for the hyper-Poisson distribution and "\code{negative.binomial}" for the negative binomial distribution. If none is specified, the default is "\code{hyper.poisson}".
#' @param method a character string, the greedy search method to be used. Possible values are "\code{hc}" for the hill-climbing and "\code{tabu}" for the tabu search. If none is specified, the default is "\code{hc}". Now, the method "\code{tabu}" is not available for the ghpd "\code{negative.binomial}".
#' @param tabu a positive integer number, the length of the tabu list for the "\code{tabu}" search. 
#' @param max.tabu a positive integer number, the iterations that the "\code{tabu}" search can perform without improving the best network score (BIC). 
#' @param maxiter an integer, the maximum number of iterations.
#' @param tol a numeric value, the tolerance for the convergence of the network score (BIC).
#' @param optim.control a list of control parameters passed to \code{optim}.
#' @param verbose a boolean value. If \code{TRUE}, progress of the algorithm is printed; otherwise the function is completely silent. 
#'
#' @return An object of class \code{linear.zidag}, a list containing the following components:\itemize{
#' \item\code{est}: a list of model parameter estimates for the linear ZiG-DAG, of which the component "\code{E}" gives the adjacency matrix of the estimated DAG for the linear ZiG-DAG.\cr
#' \item\code{bic}: the Bayesian Information Criterion for the estimated linear ZiG-DAG model.
#' \item\code{iter}: the number of iterations of the score-based greedy search algorithm used. 
#' }
#' 
#' @export
#'
#' @examples
linear.zigdag = function(dat, start = NULL, ghpd = "hyper.poisson", method = "hc", tabu = 10, max.tabu = tabu, maxiter = 500, tol = .Machine$double.eps^0.25, optim.control = list(), verbose = FALSE)
{
   if (is.null(start))
   {
      start = matrix(0, ncol(dat), ncol(dat))
   }
   
   if (ghpd == "hyper.poisson")
   {
      if (method == "hc")
      {
         out = hc_linear_ZIHP(dat = dat, starting.dag = start, maxiter = maxiter, tol = tol, optim.control = optim.control, verbose = verbose)
      } else if (method == "tabu")
      {
         out = ts_linear_ZIHP(dat = dat, starting.dag = start, tabu = tabu, max.tabu = max.tabu, maxiter = maxiter, tol = tol, optim.control = optim.control, verbose = verbose)
      }
   } else if (ghpd == "negative.binomial")
   {
      out = hc_linear_ZINB(dat = dat, starting.dag = start, maxiter = maxiter, tol = tol, optim.control = list(), verbose = verbose)
   }
   
   class(out) = "linear.zigdag"
   return(out)
}


# hill climbing for linear ZIHPBN
hc_linear_ZIHP = function(dat, starting.dag, maxiter = 500, tol = .Machine$double.eps^0.25, optim.control = list(), verbose = FALSE)
{
   n = nrow(dat) 
   p = ncol(dat)
   
   control = list(fnscale = -1, maxit = 10000, reltol = 1.0e-8)
   control[names(optim.control)] = optim.control
   
   # fit the linear ZIHPBN given the starting DAG
   bic_curr = rep(NA, p)
   est_curr = list()
   est_curr$E      = starting.dag
   est_curr$alpha  = matrix(0, p, p)
   est_curr$beta   = matrix(0, p, p)
   est_curr$delta  = rep(0, p)
   est_curr$gamma  = rep(0, p)
   est_curr$lambda = rep(1, p)
   
   for (j in 1 : p)
   {
      pa_j    = (est_curr$E[j, ] == 1) 
      start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$lambda[j]))
      out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                      gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
      p_j = sum(pa_j)
      bic_curr[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
      
      if (p_j > 0)
      {
         est_curr$alpha[j, pa_j] = out_j$par[1 : p_j]
         est_curr$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
      }
      est_curr$delta[j]  = out_j$par[p_j + 1]
      est_curr$gamma[j]  = out_j$par[2 * p_j + 2]
      est_curr$lambda[j] = exp(out_j$par[2 * p_j + 3])
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
               start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$lambda[j]))
               out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                               gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
               
               p_j = sum(pa_j)
               bic_cand = bic_curr
               bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)

               if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
               {
                  IMPROV  = TRUE
                  bic_iter = bic_cand
                  est_iter = est_curr
                  est_iter$E = E_cand
                  est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                  est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
                  est_iter$delta[j]  = out_j$par[p_j + 1]
                  est_iter$gamma[j]  = out_j$par[2 * p_j + 2]
                  est_iter$lambda[j] = exp(out_j$par[2 * p_j + 3])
               }
            } else
            {
               E_cand = est_curr$E
               E_cand[j, k] = 0
               
               pa_j    = (E_cand[j, ] == 1) 
               start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$lambda[j]))
               out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                               gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
               
               p_j = sum(pa_j)
               bic_cand = bic_curr
               bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)

               if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
               {
                  IMPROV  = TRUE
                  bic_iter = bic_cand
                  est_iter = est_curr
                  est_iter$E = E_cand
                  est_iter$alpha[j, k] = est_iter$beta[j, k] = 0
                  if (p_j > 0)
                  {
                     est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                     est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
                  }
                  est_iter$delta[j]  = out_j$par[p_j + 1]
                  est_iter$gamma[j]  = out_j$par[2 * p_j + 2]
                  est_iter$lambda[j] = exp(out_j$par[2 * p_j + 3])
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
            start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$lambda[j]))
            out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                            gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
            pa_k    = (E_cand[k, ] == 1) 
            start_k = c(est_curr$alpha[k, pa_k], est_curr$delta[k], est_curr$beta[k, pa_k], est_curr$gamma[k], log(est_curr$lambda[k]))
            out_k   = optim(par = start_k, fn = function(z) dZIHP_cpp(z, y = dat[, k], X = dat[ , pa_k, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                            gr = function(z) gradZIHP_cpp(z, y = dat[, k], X = dat[ , pa_k, drop = FALSE]), method = "BFGS", control = control)
            
            p_j = sum(pa_j)
            p_k = sum(pa_k)
            bic_cand = bic_curr
            bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
            bic_cand[k] = -2 * out_k$value + (2 * p_k + 3) * log(n)

            if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
            {
               IMPROV  = TRUE
               bic_iter = bic_cand
               est_iter = est_curr
               est_iter$E = E_cand
               est_iter$alpha[j, k] = est_iter$beta[j, k] = 0
               if (p_j > 0)
               {
                  est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                  est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
               }
               est_iter$delta[j]  = out_j$par[p_j + 1]
               est_iter$gamma[j]  = out_j$par[2 * p_j + 2]
               est_iter$lambda[j] = exp(out_j$par[2 * p_j + 3])
               est_iter$alpha[k, pa_k] = out_k$par[1 : p_k]
               est_iter$beta[k, pa_k]  = out_k$par[(p_k + 2) : (2 * p_k + 1)]
               est_iter$delta[k]  = out_k$par[p_k + 1]
               est_iter$gamma[k]  = out_k$par[2 * p_k + 2]
               est_iter$lambda[k] = exp(out_k$par[2 * p_k + 3])
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
   return(out)
}

# hill climbing for linear ZINBBN
hc_linear_ZINB = function(dat, starting.dag, maxiter = 500, tol = .Machine$double.eps^0.25, optim.control = list(), verbose = FALSE)
{
   n = nrow(dat) 
   p = ncol(dat)
   
   control = list(fnscale = -1, maxit = 10000, reltol = 1.0e-8)
   control[names(optim.control)] = optim.control
   
   # fit the linear ZINBBN given the starting DAG
   bic_curr = rep(NA, p)
   est_curr = list()
   est_curr$E      = starting.dag
   est_curr$alpha  = matrix(0, p, p)
   est_curr$beta   = matrix(0, p, p)
   est_curr$delta  = rep(0, p)
   est_curr$gamma  = rep(0, p)
   est_curr$k      = rep(1, p)
   
   for (j in 1 : p)
   {
      pa_j    = (est_curr$E[j, ] == 1) 
      start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$k[j]))
      out_j   = optim(par = start_j, fn = function(z) dZINB_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                      gr = function(z) gradZINB_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
      p_j = sum(pa_j)
      bic_curr[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
      
      if (p_j > 0)
      {
         est_curr$alpha[j, pa_j] = out_j$par[1 : p_j]
         est_curr$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
      }
      est_curr$delta[j] = out_j$par[p_j + 1]
      est_curr$gamma[j] = out_j$par[2 * p_j + 2]
      est_curr$k[j]     = exp(out_j$par[2 * p_j + 3])
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
               start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$k[j]))
               out_j   = optim(par = start_j, fn = function(z) dZINB_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                               gr = function(z) gradZINB_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
               
               p_j = sum(pa_j)
               bic_cand = bic_curr
               bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
               
               if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
               {
                  IMPROV  = TRUE
                  bic_iter = bic_cand
                  est_iter = est_curr
                  est_iter$E = E_cand
                  est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                  est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
                  est_iter$delta[j] = out_j$par[p_j + 1]
                  est_iter$gamma[j] = out_j$par[2 * p_j + 2]
                  est_iter$k[j]     = exp(out_j$par[2 * p_j + 3])
               }
            } else
            {
               E_cand = est_curr$E
               E_cand[j, k] = 0
               
               pa_j    = (E_cand[j, ] == 1) 
               start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$k[j]))
               out_j   = optim(par = start_j, fn = function(z) dZINB_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                               gr = function(z) gradZINB_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
               
               p_j = sum(pa_j)
               bic_cand = bic_curr
               bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
               
               if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
               {
                  IMPROV  = TRUE
                  bic_iter = bic_cand
                  est_iter = est_curr
                  est_iter$E = E_cand
                  est_iter$alpha[j, k] = est_iter$beta[j, k] = 0
                  if (p_j > 0)
                  {
                     est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                     est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
                  }
                  est_iter$delta[j] = out_j$par[p_j + 1]
                  est_iter$gamma[j] = out_j$par[2 * p_j + 2]
                  est_iter$k[j]     = exp(out_j$par[2 * p_j + 3])
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
            start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$k[j]))
            out_j   = optim(par = start_j, fn = function(z) dZINB_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                            gr = function(z) gradZINB_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
            pa_k    = (E_cand[k, ] == 1) 
            start_k = c(est_curr$alpha[k, pa_k], est_curr$delta[k], est_curr$beta[k, pa_k], est_curr$gamma[k], log(est_curr$k[k]))
            out_k   = optim(par = start_k, fn = function(z) dZINB_cpp(z, y = dat[, k], X = dat[ , pa_k, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                            gr = function(z) gradZINB_cpp(z, y = dat[, k], X = dat[ , pa_k, drop = FALSE]), method = "BFGS", control = control)
            
            p_j = sum(pa_j)
            p_k = sum(pa_k)
            bic_cand = bic_curr
            bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
            bic_cand[k] = -2 * out_k$value + (2 * p_k + 3) * log(n)
            
            if (is.finite(sum(bic_cand)) & (sum(bic_cand) < sum(bic_iter) - tol))
            {
               IMPROV  = TRUE
               bic_iter = bic_cand
               est_iter = est_curr
               est_iter$E = E_cand
               est_iter$alpha[j, k] = est_iter$beta[j, k] = 0
               if (p_j > 0)
               {
                  est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                  est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
               }
               est_iter$delta[j] = out_j$par[p_j + 1]
               est_iter$gamma[j] = out_j$par[2 * p_j + 2]
               est_iter$k[j]     = exp(out_j$par[2 * p_j + 3])
               est_iter$alpha[k, pa_k] = out_k$par[1 : p_k]
               est_iter$beta[k, pa_k]  = out_k$par[(p_k + 2) : (2 * p_k + 1)]
               est_iter$delta[k] = out_k$par[p_k + 1]
               est_iter$gamma[k] = out_k$par[2 * p_k + 2]
               est_iter$k[k]     = exp(out_k$par[2 * p_k + 3])
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
   return(out)
}

# tabu search algorithm for linear ZIHPBN
ts_linear_ZIHP = function(dat, starting.dag, tabu = 10, max.tabu = tabu, maxiter = 500, tol = .Machine$double.eps^0.25, optim.control = list(), verbose = FALSE)
{
   n = nrow(dat) 
   p = ncol(dat)
   
   control = list(fnscale = -1, maxit = 10000, reltol = 1.0e-8)
   control[names(optim.control)] = optim.control
   
   # fit the linear ZIHPBN given the starting DAG
   bic_curr = rep(NA, p)
   est_curr = list()
   est_curr$E      = starting.dag
   est_curr$alpha  = matrix(0, p, p)
   est_curr$beta   = matrix(0, p, p)
   est_curr$delta  = rep(0, p)
   est_curr$gamma  = rep(0, p)
   est_curr$lambda = rep(1, p)
   
   for (j in 1 : p)
   {
      pa_j    = (est_curr$E[j, ] == 1) 
      start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$lambda[j]))
      out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                      gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
      p_j = sum(pa_j)
      
      bic_curr[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
      
      if (p_j > 0)
      {
         est_curr$alpha[j, pa_j] = out_j$par[1 : p_j]
         est_curr$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
      }
      est_curr$delta[j]  = out_j$par[p_j + 1]
      est_curr$gamma[j]  = out_j$par[2 * p_j + 2]
      est_curr$lambda[j] = exp(out_j$par[2 * p_j + 3])
   }
   
   # start the tabu search
   bic_best = bic_curr
   est_best = est_curr
   LastImprovement = 0
   list_tabu = matrix(NA, tabu, 3)
   
   for (iter in 1 : maxiter)
   {
      bic_iter  = NULL
      est_iter  = NULL
      tabu_iter = NULL
      
      # search all DAGs reachable from the current DAG
      # addition or deletion of an edge
      for (j in 1 : p)
      {
         for (k in 1 : p)
         {
            if (j == k) next
            
            if (est_curr$E[j, k] == 0)
            {
               if (!legal_tabu(j, k, 1, list_tabu)) next
               
               E_cand = est_curr$E
               E_cand[j , k] = 1
               G_cand = igraph::graph_from_adjacency_matrix(E_cand)
               if (!igraph::is.dag(G_cand)) next
               
               pa_j    = (E_cand[j, ] == 1) 
               start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$lambda[j]))
               out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                               gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
               
               p_j = sum(pa_j)
               bic_cand = bic_curr
               bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
               
               if (is.finite(sum(bic_cand)) & (is.null(bic_iter) | (sum(bic_cand) < sum(bic_iter) - tol)))
               {
                  bic_iter = bic_cand
                  
                  est_iter = est_curr
                  est_iter$E = E_cand
                  est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                  est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
                  est_iter$delta[j]  = out_j$par[p_j + 1]
                  est_iter$gamma[j]  = out_j$par[2 * p_j + 2]
                  est_iter$lambda[j] = exp(out_j$par[2 * p_j + 3])
                  
                  tabu_iter = c(j, k, 1)
               }
            } else
            {
               if (!legal_tabu(j, k, -1, list_tabu)) next
               
               E_cand = est_curr$E
               E_cand[j, k] = 0
               
               pa_j    = (E_cand[j, ] == 1) 
               start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$lambda[j]))
               out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                               gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
               
               p_j = sum(pa_j)
               bic_cand = bic_curr
               bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
               
               if (is.finite(sum(bic_cand)) & (is.null(bic_iter) | (sum(bic_cand) < sum(bic_iter) - tol)))
               {
                  bic_iter = bic_cand
                  
                  est_iter = est_curr
                  est_iter$E = E_cand
                  est_iter$alpha[j, k] = est_iter$beta[j, k] = 0
                  if (p_j > 0)
                  {
                     est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                     est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
                  }
                  est_iter$delta[j]  = out_j$par[p_j + 1]
                  est_iter$gamma[j]  = out_j$par[2 * p_j + 2]
                  est_iter$lambda[j] = exp(out_j$par[2 * p_j + 3])
                  
                  tabu_iter = c(j, k, -1)
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
            if (!legal_tabu(j, k, 0, list_tabu)) next
            
            E_cand = est_curr$E
            E_cand[j , k] = 0
            E_cand[k , j] = 1
            G_cand = igraph::graph_from_adjacency_matrix(E_cand)
            if (!igraph::is.dag(G_cand)) next
            
            pa_j    = (E_cand[j, ] == 1) 
            start_j = c(est_curr$alpha[j, pa_j], est_curr$delta[j], est_curr$beta[j, pa_j], est_curr$gamma[j], log(est_curr$lambda[j]))
            out_j   = optim(par = start_j, fn = function(z) dZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                            gr = function(z) gradZIHP_cpp(z, y = dat[, j], X = dat[ , pa_j, drop = FALSE]), method = "BFGS", control = control)
            pa_k    = (E_cand[k, ] == 1) 
            start_k = c(est_curr$alpha[k, pa_k], est_curr$delta[k], est_curr$beta[k, pa_k], est_curr$gamma[k], log(est_curr$lambda[k]))
            out_k   = optim(par = start_k, fn = function(z) dZIHP_cpp(z, y = dat[, k], X = dat[ , pa_k, drop = FALSE], lower = -.Machine$double.xmax, upper = .Machine$double.xmax), 
                            gr = function(z) gradZIHP_cpp(z, y = dat[, k], X = dat[ , pa_k, drop = FALSE]), method = "BFGS", control = control)
            
            p_j  = sum(pa_j)
            p_k  = sum(pa_k)
            bic_cand = bic_curr
            bic_cand[j] = -2 * out_j$value + (2 * p_j + 3) * log(n)
            bic_cand[k] = -2 * out_k$value + (2 * p_k + 3) * log(n)
            
            if (is.finite(sum(bic_cand)) & (is.null(bic_iter) | (sum(bic_cand) < sum(bic_iter) - tol)))
            {
               bic_iter = bic_cand
               
               est_iter = est_curr
               est_iter$E = E_cand
               est_iter$alpha[j, k] = est_iter$beta[j, k] = 0
               if (p_j > 0)
               {
                  est_iter$alpha[j, pa_j] = out_j$par[1 : p_j]
                  est_iter$beta[j, pa_j]  = out_j$par[(p_j + 2) : (2 * p_j + 1)]
               }
               est_iter$delta[j]  = out_j$par[p_j + 1]
               est_iter$gamma[j]  = out_j$par[2 * p_j + 2]
               est_iter$lambda[j] = exp(out_j$par[2 * p_j + 3])
               est_iter$alpha[k, pa_k] = out_k$par[1 : p_k]
               est_iter$beta[k, pa_k]  = out_k$par[(p_k + 2) : (2 * p_k + 1)]
               est_iter$delta[k]  = out_k$par[p_k + 1]
               est_iter$gamma[k]  = out_k$par[2 * p_k + 2]
               est_iter$lambda[k] = exp(out_k$par[2 * p_k + 3])
               
               tabu_iter = c(k, j, 0)
            }
         }
      }
      
      if (is.null(bic_iter)) break
      
      if (verbose) 
         cat("iter =", iter, "; BIC =", round(sum(bic_iter), 4), "\n")
      
      bic_curr = bic_iter
      est_curr = est_iter
      
      if (tabu > 0)
      {
         if (tabu > 1) list_tabu[2 : tabu, ] = list_tabu[1 : (tabu - 1), ]
         list_tabu[1, ] = tabu_iter
      }
      
      if (sum(bic_curr) < sum(bic_best) - tol)
      {
         bic_best = bic_curr
         est_best = est_curr
         LastImprovement = 0
      } else
      {
         LastImprovement = LastImprovement + 1
      }
      
      if (LastImprovement >= max.tabu) break
   }
   
   out = list()
   out$est  = est_best
   out$bic  = sum(bic_best)
   out$iter = iter
   return(out)
}

# check if a new operator is legal given the current tabu list
legal_tabu = function(j, k, opr, list_tabu)        # opr = 1:addition, = -1:deletion, = 0:reversal
{
   out = TRUE
   id_tabu = which(list_tabu[ , 3] == -opr)
   tabus   = list_tabu[id_tabu, -3, drop = FALSE]
   if (any(tabus[ , 1] == j & tabus[ , 2] == k)) out = FALSE
   
   return(out)
}
