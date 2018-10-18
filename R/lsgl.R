#' @title Linear Multiple Output Using Sparse Group Lasso.
#'
#' @description Simultaneous feature selection and parameter estimation for linear multiple output.
#' The algorithm finds the sparse group lasso penalized maximum likelihood estimator.
#' This result in feature and parameter selection, and parameter estimation.
#' Use of parallel computing for cross validation and subsampling is supported through the 'foreach' and 'doParallel' packages.
# 'Development version is on GitHub, please report package issues on GitHub.
#'
#' @details
#' The sparse gorup lasso linear mltiple output estimator is defined as
#' \deqn{\frac{1}{N}\|Y-X\beta\|_F^2 + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
#' where \eqn{\|\cdot\|_F} is the frobenius norm.
#' The vector \eqn{\beta^{(J)}} denotes the parameters associated with the \eqn{J}'th group of features.
#' The group weights are denoted by \eqn{\gamma \in [0,\infty)^m} and the parameter weights by \eqn{\xi \in [0,\infty)^n}.
#'
#' @author Martin Vincent \email{martin.vincent.dk@gmail.com}
#'
#' @examples
#' # NOTE
#' @docType package
#' @name lsgl-package
#' @importFrom tools assertWarning
#' @useDynLib lsgl, .registration=TRUE
NULL

#' Airline Ticket Prices.
#'
#' @format A design matrix and a response matrix
#' \describe{
#'   \item{X}{design matrix}
#'   \item{Y}{response matrix}
#' }
#' @name AirlineTicketPrices
#' @docType data
#' @keywords data
NULL

#' Design matrix
#' @name X
#' @docType data
#' @keywords data
NULL

#' Response matrix
#' @name Y
#' @docType data
#' @keywords data
NULL
