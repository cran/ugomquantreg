#' @name ammonia
#' @aliases ammonia
#' 
#' @title Ammonia oxidized to acid nitric data set
#'
#' @description The data come from experiments with a plant where ammonia is oxidized to acid nitric.
#'
#' @format A data-frame with 21 observations and 4 columns:
#'
#' \itemize{
#' \item \code{stackloss}: the percentage of ammonia lost.
#' \item \code{airflow}: the air flow to the plant.
#' \item \code{watertemp}: the cooling water inlet temperature.
#' \item \code{acidconc}: the acid concentration.
#' }
#' 
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' 
#' @author Bruna Alves \email{pg402900@uem.br}
#' 
#' @usage data(ammonia, package = "ugomquantreg")
#'
#' @references 
#' 
#' Brownlee, K. A., (1965). Statistical Theory and Methodology in Science and Engineering. \emph{New York: John Wiley & Sons}.
#' 
#' Yu, K., and Moyeed, R. A., (2001). Bayesian quantile regression. \emph{Statistics and Probability Letters}, \bold{54}(4) 437--447.
#' 
#' @source \url{https://support.sas.com/rnd/app/stat/examples/BayesQuantile/quantile.htm}
#' 
#' @examples 
#' data(ammonia, package = "ugomquantreg")
#' 
#' library(gamlss)
#' 
#' tau <- 0.50
#' fit.logit <- gamlss(stackloss ~ airflow + watertemp + acidconc, data = ammonia, 
#' family = UGOM(sigma.link="identity"))
#' 
#' tau <- 0.50
#' fit.probit <- gamlss(stackloss ~ airflow + watertemp + acidconc, 
#' data = ammonia, family = UGOM(mu.link = "probit", sigma.link = "log"))
#' 
#' fittaus <- lapply(c(0.10, 0.25, 0.50, 0.75, 0.90), function(Tau){
#'  tau <<- Tau;
#'  gamlss(stackloss ~ airflow + watertemp + acidconc, data = ammonia, 
#'  family = UGOM(mu.link = "logit", sigma.link = "log"))
#' })
#' 
#' sapply(fittaus, coef)
"ammonia"

