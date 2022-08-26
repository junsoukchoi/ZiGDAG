#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

double F11_cpp(double z, double a, double b, int maxiter, double tol)
{
  double fac = 1;
  double temp = fac;
  double series = 0;
  
  for (int n = 1; n < maxiter + 1; n++)
  {
    fac = fac * (a / b) * (z / n);
    if (Rcpp::traits::is_nan<REALSXP>(fac)) fac = 0;
    series = temp + fac;
    
    if (Rcpp::traits::is_infinite<REALSXP>(series) | (std::abs(series - temp) < tol))
    {
      return(series);
    }

    temp = series;
    a++;
    b++;
  }
  
  return(series);
}

double gradF11b_cpp(double z, double  b, int maxiter, double tol)
{
  double con  = R::digamma(b) * F11_cpp(z, 1, b, maxiter, tol);
  double fac  = R::digamma(b);
  double temp = fac;
  double series = 0;
  
  for (int n = 1; n < maxiter + 1; n++)
  {
    fac = fac * R::digamma(b + 1) / R::digamma(b) * z / b;
    if (Rcpp::traits::is_nan<REALSXP>(fac)) fac = 0;
    series = temp + fac;
    
    if (Rcpp::traits::is_infinite<REALSXP>(series) | (std::abs(series - temp) < tol))
    {
      return(con - series);
    }
    
    temp = series;
    b++;
  }
  
  return(con - series);
}

arma::vec F11_for_each(const arma::vec & z, double a, double b, int maxiter, double tol)
{
  int n = z.n_elem;
  arma::vec val(n);
  
  for (int i = 0; i < n; i++)
  {
    val(i) = F11_cpp(z(i), a, b, maxiter, tol);
  }
  
  return(val);
}

arma::vec gradF11b_for_each(const arma::vec & z, double b, int maxiter, double tol)
{
  int n = z.n_elem;
  arma::vec val(n);
  
  for (int i = 0; i < n; i++)
  {
    val(i) = gradF11b_cpp(z(i), b, maxiter, tol);
  }
     
  return(val);
}

arma::vec log_ascfacto(double z, const arma::vec & n)
{
  int p = n.n_elem;
  arma::vec out(p, arma::fill::zeros);
  
  for (int i = 0; i < p; i++)
  {
    for (int j = 0; j < n(i); j++)
    {
      out(i) += log(z + j);
    }
  }
  
  return(out);
}

arma::vec digamma_arma(const arma::vec & z)
{
  Rcpp::NumericVector v = Rcpp::digamma(Rcpp::NumericVector(z.begin(), z.end()));
  arma::vec out(v.begin(), v.size(), false);
  return(out);
}

// [[Rcpp::export]]
double dZIHP_cpp(const arma::vec & param, const arma::vec & y, const arma::mat & X, double lower, double upper)
{
  double n = X.n_rows;
  double p = X.n_cols;
  arma::mat X1 = arma::join_rows(X, arma::ones<arma::vec>(n));
  arma::vec Pi(n);
  Pi = arma::exp(X1 * param.subvec(0, p));
  Pi = Pi / (1 + Pi);
  Pi.elem(arma::find_nonfinite(Pi)).fill(1);
  arma::vec Mu(n);
  Mu = arma::exp(X1 * param.subvec(p + 1, 2 * p + 1));
  double Lambda = std::exp(param(2 * p + 2));
  
  arma::uvec Y0 = arma::find(abs(y) <  1.0e-8);
  arma::uvec Y1 = arma::find(abs(y) >= 1.0e-8);
  arma::vec  eval_F11 = F11_for_each(Mu, 1, Lambda, 10000, 1.0e-8);

  double llik = arma::accu(arma::log(Pi.elem(Y0) + (1 - Pi.elem(Y0)) / eval_F11.elem(Y0))) + 
    arma::accu(arma::log(1 - Pi.elem(Y1)) - log_ascfacto(Lambda, y.elem(Y1)) - arma::log(eval_F11.elem(Y1)) + y.elem(Y1) % arma::log(Mu.elem(Y1)));
  
  if (Rcpp::traits::is_infinite<REALSXP>(llik) & (llik < 0))
  {
    return(lower);
  } else if (Rcpp::traits::is_infinite<REALSXP>(llik) & (llik > 0))
  {
    return(upper);
  } else
  {
    return(llik);
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector gradZIHP_cpp(const arma::vec & param, const arma::vec & y, const arma::mat & X)
{
  double n = X.n_rows;
  double p = X.n_cols;
  arma::mat X1 = arma::join_rows(X, arma::ones<arma::vec>(n));
  arma::vec Pi(n);
  Pi = arma::exp(X1 * param.subvec(0, p));
  Pi = Pi / (1 + Pi);
  Pi.elem(arma::find_nonfinite(Pi)).fill(1);
  arma::vec Mu(n);
  Mu = arma::exp(X1 * param.subvec(p + 1, 2 * p + 1));
  double Lambda = std::exp(param(2 * p + 2));
  
  arma::uvec Y0 = arma::find(abs(y) <  1.0e-8);
  arma::uvec Y1 = arma::find(abs(y) >= 1.0e-8);
  arma::mat  X10 = X1.rows(Y0);
  arma::mat  X11 = X1.rows(Y1);
  arma::vec  eval1_F11 = F11_for_each(Mu, 1, Lambda, 10000, 1.0e-8);
  arma::vec  eval2_F11 = F11_for_each(Mu, 2, Lambda + 1, 10000, 1.0e-8);
  arma::vec  eval_gradF11b = gradF11b_for_each(Mu, Lambda, 10000, 1.0e-8);
  arma::vec  dZIHP0 = Pi.elem(Y0) + (1 - Pi.elem(Y0)) / eval1_F11.elem(Y0);
  
  arma::rowvec Grad(2 * p + 3);
  Grad.subvec(0, p) = arma::sum(X10.each_col() % ((1 - Pi.elem(Y0)) % (Pi.elem(Y0) - Pi.elem(Y0) / eval1_F11.elem(Y0)) / dZIHP0), 0) - arma::sum(X11.each_col() % Pi.elem(Y1), 0);
  Grad.subvec(p + 1, 2 * p + 1) = -arma::sum(X10.each_col() % (Mu.elem(Y0) % (1 - Pi.elem(Y0)) / Lambda % eval2_F11.elem(Y0) / arma::square(eval1_F11.elem(Y0)) / dZIHP0), 0) -
    arma::sum(X11.each_col() % (Mu.elem(Y1) / Lambda % eval2_F11.elem(Y1) / eval1_F11.elem(Y1) - y.elem(Y1)), 0);
  Grad(2 * p + 2) = (-arma::accu((1 - Pi.elem(Y0)) % eval_gradF11b.elem(Y0) / arma::square(eval1_F11.elem(Y0)) / dZIHP0) - 
    arma::accu(eval_gradF11b.elem(Y1) / eval1_F11.elem(Y1) + digamma_arma(Lambda + y.elem(Y1)) - R::digamma(Lambda))) * Lambda;
  
  Rcpp::NumericVector out(Grad.begin(), Grad.end());
  return(out);
}

// [[Rcpp::export]]
double dZINB_cpp(const arma::vec & param, const arma::vec & y, const arma::mat & X, double lower, double upper)
{
  double n = X.n_rows;
  double p = X.n_cols;
  arma::mat X1 = arma::join_rows(X, arma::ones<arma::vec>(n));
  arma::vec pi(n);
  pi = arma::exp(X1 * param.subvec(0, p));
  pi = pi / (1 + pi);
  pi.elem(arma::find_nonfinite(pi)).fill(1);
  arma::vec q(n);
  q = arma::exp(X1 * param.subvec(p + 1, 2 * p + 1));
  q = q / (1 + q);
  q.elem(arma::find_nonfinite(q)).fill(1);
  double k = std::exp(param(2 * p + 2));
  
  arma::uvec Y0 = arma::find(abs(y) <  1.0e-8);
  arma::uvec Y1 = arma::find(abs(y) >= 1.0e-8);
  
  double llik = arma::accu(arma::log(pi.elem(Y0) + (1 - pi.elem(Y0)) % arma::pow(1 - q.elem(Y0), k))) + 
    arma::accu(arma::log(1 - pi.elem(Y1)) + log_ascfacto(k, y.elem(Y1)) - arma::lgamma(y.elem(Y1) + 1) + y.elem(Y1) % arma::log(q.elem(Y1)) + k * arma::log(1 - q.elem(Y1)));
  
  if (Rcpp::traits::is_infinite<REALSXP>(llik) & (llik < 0))
  {
    return(lower);
  } else if (Rcpp::traits::is_infinite<REALSXP>(llik) & (llik > 0))
  {
    return(upper);
  } else
  {
    return(llik);
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector gradZINB_cpp(const arma::vec & param, const arma::vec & y, const arma::mat & X)
{
  double n = X.n_rows;
  double p = X.n_cols;
  arma::mat X1 = arma::join_rows(X, arma::ones<arma::vec>(n));
  arma::vec pi(n);
  pi = arma::exp(X1 * param.subvec(0, p));
  pi = pi / (1 + pi);
  pi.elem(arma::find_nonfinite(pi)).fill(1);
  arma::vec q(n);
  q = arma::exp(X1 * param.subvec(p + 1, 2 * p + 1));
  q = q / (1 + q);
  q.elem(arma::find_nonfinite(q)).fill(1);
  double k = std::exp(param(2 * p + 2));
  
  arma::uvec Y0 = arma::find(abs(y) <  1.0e-8);
  arma::uvec Y1 = arma::find(abs(y) >= 1.0e-8);
  arma::mat  X10 = X1.rows(Y0);
  arma::mat  X11 = X1.rows(Y1);
  arma::vec  dZINB0 = pi.elem(Y0) + (1 - pi.elem(Y0)) % arma::pow(1 - q.elem(Y0), k);
  
  arma::rowvec Grad(2 * p + 3);
  Grad.subvec(0, p) = arma::sum(X10.each_col() % ((1 - pi.elem(Y0)) % (pi.elem(Y0) - pi.elem(Y0) % arma::pow(1 - q.elem(Y0), k)) / dZINB0), 0) - 
    arma::sum(X11.each_col() % pi.elem(Y1), 0);
  Grad.subvec(p + 1, 2 * p + 1) = -arma::sum(X10.each_col() % ((1 - pi.elem(Y0)) % (k * arma::pow(1 - q.elem(Y0), k - 1)) % q.elem(Y0) % (1 - q.elem(Y0)) / dZINB0), 0) -
    arma::sum(X11.each_col() % ((k + y.elem(Y1)) % q.elem(Y1) - y.elem(Y1)), 0);
  Grad(2 * p + 2) = (arma::accu((1 - pi.elem(Y0)) % arma::pow(1 - q.elem(Y0), k) % arma::log(1 - q.elem(Y0)) / dZINB0) + 
    arma::accu(digamma_arma(k + y.elem(Y1)) - R::digamma(k) + arma::log(1 - q.elem(Y1)))) * k;
  
  Rcpp::NumericVector out(Grad.begin(), Grad.end());
  return(out);
}