#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf unit-Gompertz

inline double logpdf_ugompertz(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (0.1e1 / pow(mu, theta) - 0.1e1);
  double t1 = log(alpha);
  double t2 = log(theta);
  double t3 = lnx;
  double t5 = pow(x, theta);
  return(t1 + t2 - theta * t3 + alpha - alpha / t5 - t3);
}

// [[Rcpp::export]]
NumericVector cpp_dugompertz(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const NumericVector tau,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();
  const int ntau = tau.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_ugompertz(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau[i % ntau]);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf unit-gompertz
inline double cdf_ugompertz(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (0.1e1 / pow(mu, theta) - 0.1e1);
  double t1 = pow(x, -theta);
  return(exp(alpha - alpha * t1));
}

// [[Rcpp::export]]
NumericVector cpp_pugompertz(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const NumericVector tau,
                             const bool lowertail = true,
                             const bool logprob = false)
{
  const int n = x.length(); 
  NumericVector out(n);
  const int nmu = mu.length(); 
  const int nth = theta.length(); 
  const int ntau = tau.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_ugompertz(x[i], mu[i % nmu], theta[i % nth], tau[i % ntau]);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);
  return(out);
}

// inv-cdf unit-gompertz
inline double invcdf_ugompertz(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (0.1e1 / pow(mu, theta) - 0.1e1);
  double t1 = log(x);
  double t6 = pow((alpha - t1) / alpha, 0.1e1 / theta);
  return(0.1e1 / t6);
}

// [[Rcpp::export]]
NumericVector cpp_qugompertz(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const NumericVector tau,
                             const bool lowertail = true,
                             const bool logprob = false)
{
  const int n = x.length(); 
  NumericVector out(n);
  const int nmu = mu.length(); 
  const int nth = theta.length();
  const int ntau = tau.length();

  if(lowertail)
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ugompertz(x[i], mu[i % nmu], theta[i % nth], tau[i % ntau]);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ugompertz(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau[i % ntau]);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// cdf unit-gompertz
inline double DMU(double x, double mu, double sigma, double nu)
{
  double t2 = pow(mu, sigma);
  double t4 = pow(-0.1e1 + t2, 0.2e1);
  double t9 = pow(mu, 0.1e1 + sigma);
  double t13 = pow(mu, 0.1e1 + 0.2e1 * sigma);
  double t14 = log(nu);
  double t16 = pow(x, -sigma);
  return(-0.1e1 / mu / t4 * (-sigma + sigma * t2 - mu + 0.2e1 * t9 - t13 - t14 * t2 * t16 * sigma));
}

// [[Rcpp::export]]
NumericVector cpp_dmu(const NumericVector x,
                      const NumericVector mu,
                      const NumericVector sigma,
                      const NumericVector nu)
{
  const int n = x.length(); 
  NumericVector out(n);
  for(int i = 0; i < n; i++)
    out[i] = DMU(x[i], mu[i], sigma[i], nu[i]);
  
  return(out);
}
