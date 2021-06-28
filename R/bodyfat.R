#' @name bodyfat
#' @aliases bodyfat
#' 
#' @title Percentage of body fat data set
#'
#' @description The body fat percentage of individuals assisted in a public hospital in Curitiba, Paraná, Brazil.
#'
#' @format A data-frame with 298 observations and 9 columns:
#'
#' \itemize{
#' \item \code{ARMS}: arms fat percentage.
#' \item \code{LEGS}: legs fat percentage.
#' \item \code{BODY}: body fat percentage.
#' \item \code{ANDROID}: android fat percentage.
#' \item \code{GYNECOID}: ginecoid fat percentage.
#' \item \code{AGE}: age of individuals.
#' \item \code{BMI}: body mass index.
#' \item \code{SEX}: 1 for female, 2 for male.
#' \item \code{IPAQ}: 0 for IPAQ = sedentary, 1 for IPAQ = insufficiently active and 2 for IPAQ = active.
#' }
#' 
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' 
#' @author Bruna Alves \email{pg402900@uem.br}
#' 
#' @usage data(bodyfat, package = "ugomquantreg")
#'
#' @references 
#' 
#' Mazucheli, J., Leiva, V., Alves, B., and Menezes A. F. B., (2021). A new quantile regression for modeling bounded data under a unit Birnbaum-Saunders distribution with applications in medicine and politics. \emph{Symmetry}, \bold{13}(4) 1--21.
#' 
#' Petterle, R. R., Bonat, W. H., Scarpin, C. T., Jonasson, T., and Borba, V. Z. C., (2020). Multivariate quasi-beta regression models for continuous bounded data. \emph{The International Journal of Biostatistics}, 1--15, (preprint).
#' 
#' @source \url{http://www.leg.ufpr.br/doku.php/publications:papercompanions:multquasibeta}
#' 
#' @examples 
#' data(bodyfat, package = "ugomquantreg")
#' 
#' library(gamlss)
#' 
#' tau <- 0.50
#' fit.logit <- gamlss(ARMS ~ AGE + I(BMI / 100) + as.factor(SEX) + as.factor(IPAQ), 
#' data = bodyfat, family = UGOM(mu.link = "logit", sigma.link = "log"))
#' 
#' tau <- 0.50;
#' fit.probit <- gamlss(ARMS ~ AGE + I(BMI / 100) + as.factor(SEX) + as.factor(IPAQ), 
#' data = bodyfat, family = UGOM(mu.link = "probit", sigma.link = "log"))
#' 
#'
"bodyfat"

