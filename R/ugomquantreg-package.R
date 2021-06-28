#' @docType package
#' @name ugomquantreg-package
#' @aliases ugomquantreg-package
#'
#' @title Overview of the ugomquantreg package
#'
#' @description The \pkg{ugomquantreg} package implements the probability density function, quantile function, cumulative distribution function and random number generation function for unit-Gompertz distribution parameterized as a function of its \eqn{\tau}-th quantile, \eqn{0 <\tau<1}. Some function are written in \proglang{C++} using \pkg{Rcpp}. 
#' 
#' @details 
#' 
#' \code{\link[ugomquantreg]{ammonia}}: Ammonia oxidized to acid nitric data set.
#' 
#' \code{\link[ugomquantreg]{bodyfat}}: Body fat data set.
#' 
#' \code{\link[ugomquantreg]{UGOM}}: For quantile modeling (con/in)ditional on covariate(s).
#' 
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' 
#' @author Bruna Alves \email{pg402900@uem.br}
#' 
#'
#' @useDynLib ugomquantreg
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
    library.dynam.unload("ugomquantreg", libpath)
}
