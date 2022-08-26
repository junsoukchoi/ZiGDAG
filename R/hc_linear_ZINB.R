# hill climbing algorithm for linear ZINBBN
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
