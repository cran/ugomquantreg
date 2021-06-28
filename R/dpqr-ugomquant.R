#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist checklink
#' @importFrom pracma incgam
#' @importFrom stats dnorm qnorm runif
#' @name UGOM
#' @aliases UGOM dUGOM pUGOM qUGOM rUGOM
#'
#' @title The unit-Gompertz distribution - quantile parameterization
#'
#' @description The function \code{UGOM()} define the unit-Gompertz distribution for a \code{gamlss.family} object to be used in GAMLSS fitting. \code{UGOM()} has the \eqn{\tau}-th quantile equal to the parameter mu and sigma as the shape parameter. The functions \code{dUGOM}, \code{pUGOM}, \code{qUGOM} and \code{rUGOM} define the density, distribution function, quantile function and random generation for unit-Gompertz distribution.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' 
#' @author Bruna Alves \email{pg402900@uem.br}
#' 
#' @references
#'  
#' Hastie, T. J. and Tibshirani, R. J. (1990). \emph{Generalized Additive Models}. Chapman and Hall, London.
#' 
#' Mazucheli, J., Alve, B. (2021). The Unit-Gompertz quantile regression model for bounded responses. \emph{preprint}, \bold{0}(0), 1-20.
#' 
#' Mazucheli, J., Menezes, A. F. and Dey S. (2019). Unit-Gompertz distribution with applications. \emph{Statistica}, \bold{79}(1), 25--43.
#' 
#' Rigby, R. A. and  Stasinopoulos, D. M. (2005). Generalized additive models for location, scale and shape (with discussion). \emph{Applied. Statistics}, \bold{54}(3), 507--554.
#' 
#' Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z. and De Bastiani, F. (2019). \emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}. Chapman and Hall/CRC.
#' 
#' Stasinopoulos, D. M. and Rigby, R. A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R. \emph{Journal of Statistical Software}, \bold{23}(7), 1--45.
#' 
#' Stasinopoulos, D. M., Rigby, R. A., Heller, G., Voudouris, V. and De Bastiani F. (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC.  
#'  
#' @param x,q vector of quantiles on the (0,1) interval.
#' @param p vector of probabilities.
#' @param n the number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq{x})} are returned, otherwise \eqn{P(X > x)}.
#' @param mu.link the mu link function with default logit.
#' @param sigma.link the sigma link function with default logit.
#' @param mu vector of quantile parameter values.
#' @param sigma vector of shape parameter values.
#' @param tau the \eqn{\tau}-th fixed quantile in \[d-p-q-r\]-UGOM function.
#'
#' @return \code{UGOM()} return a gamlss.family object which can be used to fit a unit-Gompertz distribution by gamlss() function.
#' 
#' @note Note that for \code{UGOM()}, mu is the \eqn{\tau}-th quantile and sigma a shape parameter. The \code{\link[gamlss]{gamlss}} function is used for parameters estimation.
#' 
#' @details
#' Probability density function 
#' \deqn{f\left( {x\mid \mu ,\sigma ,\tau } \right)=\left ( \frac{\log \left( \tau \right) }{1-\mu ^{-\sigma }} \right ) \sigma x^{-\left( 1+\sigma \right) }\exp \left[ \left ( \frac{\log \left( \tau \right) }{1-\mu ^{-\sigma }} \right ) \left( 1-x^{-\sigma }\right) \right] }
#' 
#' Cumulative distribution function
#' \deqn{F\left({x\mid \mu ,\sigma ,\tau } \right) = \exp \left[ \left ( \frac{\log \left( \tau \right) }{1-\mu ^{-\sigma }} \right ) \left( 1-x^{-\sigma }\right) \right]}
#' 
#' Mean 
#' \deqn{E(X)=\left( \frac{\log \left( \tau \right) }{1-\mu ^{-\sigma }}\right) ^{\frac{1}{\theta }}\exp \left(\frac{\log \left( \tau \right) }{1-\mu ^{-\sigma }}\right)\Gamma \left( \frac{\sigma -1}{\sigma },\frac{\log \left( \tau \right) }{  1-\mu ^{-\sigma }}\right)}
#' 
#' where \eqn{0 < (x, \mu)<1}, \eqn{\mu} is, for a fixed and known value of \eqn{\tau}, the \eqn{\tau}-th quantile, \eqn{\sigma} is the shape parameter and \eqn{\Gamma(a, b)} is the upper incomplete gamma function.

#' @examples
#'
#' set.seed(123)
#' x <- rUGOM(n = 1000, mu = 0.50, sigma = 1.69, tau = 0.50)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], length.out = 1000)
#' 
#' hist(x, prob = TRUE, main = 'unit-Gompertz')
#' lines(S, dUGOM(x = S, mu = 0.50, sigma = 1.69, tau = 0.50), col = 2)
#' 
#' plot(ecdf(x))
#' lines(S, pUGOM(q = S, mu = 0.50, sigma = 1.69, tau = 0.50), col = 2)
#' 
#' plot(quantile(x, probs = S), type = "l")
#' lines(qUGOM(p = S, mu = 0.50, sigma = 1.69, tau = 0.50), col = 2)
#' 
#' library(gamlss)
#' set.seed(123)
#' data <- data.frame(y =  rUGOM(n = 100, mu = 0.5, sigma = 2.0, tau = 0.5))
#' 
#' tau <- 0.50
#' fit <- gamlss(y ~ 1, data = data, family = UGOM)
#' 
#' set.seed(123)
#' n <- 100
#' x <- rbinom(n, size = 1, prob = 0.5)
#' eta <- 0.5 + 1 * x;
#' mu <- 1 / (1 + exp(-eta));
#' sigma <- 1.5;
#' y <- rUGOM(n, mu, sigma, tau = 0.5)
#' data <- data.frame(y, x)
#' 
#' tau <- 0.50
#' fit <- gamlss(y ~ x, data = data, family = UGOM(mu.link = "logit", sigma.link = "log"))
##################################################
#' @rdname UGOM
#' @export
#
dUGOM <- function (x, mu, sigma, tau = 0.5, log = FALSE)
{
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, sigma > 0, tau > 0, tau < 1);
  cpp_dugompertz(x, mu, sigma, tau, log[1L]);
}
##################################################
#' @rdname UGOM
#' @export
#'
pUGOM <- function (q, mu, sigma, tau = 0.50, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, sigma > 0, tau > 0, tau < 1);
  cpp_pugompertz(q, mu, sigma, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname UGOM
#' @export
#'
qUGOM <- function(p, mu, sigma, tau = 0.50, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, sigma > 0, tau > 0, tau < 1);
  cpp_qugompertz(p, mu, sigma, tau, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname UGOM
#' @export
#'
rUGOM <- function(n, mu, sigma, tau = 0.50)
{
  cpp_qugompertz(runif(n), mu, sigma, tau, TRUE, FALSE)
}
##################################################
#' @rdname UGOM
#' @export
#' 
UGOM <- function (mu.link = "logit", sigma.link = "log")
{
  mstats <- checklink("mu.link", "UGOM", substitute(mu.link), c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "UGOM", substitute(sigma.link), c("log", "identity", "own"))
  structure(
    list(family     = c("UGOM", "UGOMQ"),
         parameters = list(mu = TRUE, sigma = TRUE),
         nopar      = 2,
         type       = "Continuous",
         mu.link    = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         mu.linkfun = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
         mu.linkinv = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
              mu.dr = mstats$mu.eta,
           sigma.dr = dstats$mu.eta,
               dldm = function(y, mu, sigma, tau){
                 t1 = mu ^ sigma;
                 t2 = 0.1e1 / t1;
                 t3 = t2 * sigma;
                 t4 = 0.1e1 / mu;
                 t5 = t2 - 0.1e1;
                 t9 = log(tau);
                 t10 = t5 * t5;
                 t12 = t9 / t10;
                 t13 = t3 * t4;
                 t15 = y ^ sigma;
                 return(t3 * t4 / t5 - t12 * t13 + t12 / t15 * t13);
               },
             d2ldm2 = function(y, mu, sigma, tau){
               t1 = mu ^ sigma;
               t2 = 0.1e1 / t1;
               t3 = sigma * sigma;
               t4 = t2 * t3;
               t5 = mu * mu;
               t6 = 0.1e1 / t5;
               t7 = t2 - 0.1e1;
               t9 = t6 / t7;
               t11 = t2 * sigma;
               t13 = t1 * t1;
               t15 = 0.1e1 / t13 * t3;
               t16 = t7 * t7;
               t17 = 0.1e1 / t16;
               t20 = log(tau);
               t23 = t20 / t16 / t7;
               t24 = t15 * t6;
               t27 = t20 * t17;
               t28 = t4 * t6;
               t30 = t11 * t6;
               t32 = y ^ sigma;
               t33 = 0.1e1 / t32;
               t37 = t27 * t33;
               return(-t4 * t9 - t9 * t11 + t15 * t6 * t17 - 0.2e1 * t23 * t24 + t27 * t28 + t27 * t30 + 0.2e1 * t23 * t33 * t24 - t37 * t28 - t37 * t30);
             },
               dldd = function(y, mu, sigma, tau){
                 t1 = mu ^ sigma;
                 t2 = 0.1e1 / t1;
                 t3 = log(mu);
                 t4 = t2 * t3;
                 t5 = t2 - 0.1e1;
                 t6 = 0.1e1 / t5;
                 t9 = log(y);
                 t10 = log(tau);
                 t11 = t5 * t5;
                 t13 = t10 / t11;
                 t15 = y ^ sigma;
                 t16 = 0.1e1 / t15;
                 return(t4 * t6 + 0.1e1 / sigma - t9 - t13 * t4 + t13 * t16 * t2 * t3 - t10 * t6 * t16 * t9);
               },
      d2ldmdd = function(y, mu, sigma, tau){
        t1 = mu ^ sigma;
        t2 = 0.1e1 / t1;
        t3 = log(mu);
        t4 = t2 * t3;
        t5 = t2 - 0.1e1;
        t6 = 0.1e1 / t5;
        t8 = 0.1e1 / mu;
        t11 = t2 * t8;
        t13 = t1 * t1;
        t14 = 0.1e1 / t13;
        t15 = t14 * t3;
        t16 = t5 * t5;
        t17 = 0.1e1 / t16;
        t21 = log(tau);
        t24 = t21 / t16 / t5;
        t27 = t3 * sigma * t8;
        t30 = t21 * t17;
        t34 = y ^ sigma;
        t35 = 0.1e1 / t34;
        t37 = sigma * t8;
        t41 = t30 * t35;
        t47 = log(y);
        return(-t4 * t6 * sigma * t8 + t11 * t6 + t15 * t17 * sigma * t8 - 0.2e1 * t24 * t14 * t27 + t30 * t2 * t27 - t30 * t11 + 0.2e1 * t24 * t35 * t15 * t37 - t41 * t4 * t37 + t30 * t35 * t2 * t8 - t41 * t47 * t2 * t37);
        
        },
      d2ldd2 = function(y, mu, sigma, tau){
        t1 = mu ^ sigma;
        t2 = 0.1e1 / t1;
        t3 = log(mu);
        t4 = t3 * t3;
        t5 = t2 * t4;
        t6 = t2 - 0.1e1;
        t7 = 0.1e1 / t6;
        t9 = t1 * t1;
        t10 = 0.1e1 / t9;
        t11 = t10 * t4;
        t12 = t6 * t6;
        t13 = 0.1e1 / t12;
        t15 = sigma * sigma;
        t17 = log(tau);
        t20 = t17 / t12 / t6;
        t23 = t17 * t13;
        t25 = y ^ sigma;
        t26 = 0.1e1 / t25;
        t33 = log(y);
        t41 = t33 * t33;
        return(-t5 * t7 + t11 * t13 - 0.1e1 / t15 - 0.2e1 * t20 * t11 + t23 * t5 + 0.2e1 * t20 * t26 * t10 * t4 - 0.2e1 * t23 * t26 * t2 * t3 * t33 - t23 * t26 * t2 * t4 + t17 * t7 * t26 * t41);
               },
        G.dev.incr = function(y, mu, sigma, tau, w, ...) -2 * dUGOM(y, mu, sigma, tau, log = TRUE),
             rqres = expression(rqres(pfun = "pUGOM", type = "Continuous", y = y, mu = mu, sigma = sigma, tau = tau)),
        mu.initial = expression({mu <- (y + mean(y))/2}),
     sigma.initial = expression({sigma <- rep(1.0, length(y))}),
          mu.valid = function(mu) all(mu > 0 & mu < 1),
       sigma.valid = function(sigma) all(sigma > 0),
           y.valid = function(y) all(y > 0 & y < 1),
              mean = NULL,
     variance = NULL),class = c("gamlss.family","family"))
}
