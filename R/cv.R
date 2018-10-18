#
#     Description of this R script:
#     R interface for linear multiple output sparse group lasso routines.
#
#     Intended for use with R.
#     Copyright (C) 2014 Martin Vincent
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' @title Cross Validation
#' @description
#' Linear multiple output cross validation using multiple possessors
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param y response matrix, matrix of size \eqn{N \times K}.
#' @param intercept should the model include intercept parameters.
#' @param weights sample weights, vector of size \eqn{N \times K}.
#' @param grouping grouping of features, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the feature.
#' @param groupWeights the group weights, a vector of length \eqn{m} (the number of groups).
#' @param parameterWeights a matrix of size \eqn{K \times p}.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda lambda.min relative to lambda.max or the lambda sequence for the regularization path.
#' @param d length of lambda sequence (ignored if \code{length(lambda) > 1})
#' @param fold the fold of the cross validation, an integer larger than 1 and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' @param cv.indices a list of indices of a cross validation splitting.
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param max.threads Deprecated (will be removed in 2018),
#' instead use \code{use_parallel = TRUE} and registre parallel backend (see package 'doParallel').
#' The maximal number of threads to be used.
#' @param use_parallel If \code{TRUE} the \code{foreach} loop will use \code{\%dopar\%}. The user must registre the parallel backend.
#' @param algorithm.config the algorithm configuration to be used.
#' @return
#' \item{Yhat}{the cross validation estimated response matrix}
#' \item{Y.true}{the true response matrix, this is equal to the argument \code{y}}
#' \item{cv.indices}{the cross validation splitting used}
#' \item{features}{number of features used in the models}
#' \item{parameters}{number of parameters used in the models.}
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K
#'
#' N <- 50 #number of samples
#' p <- 25 #number of features
#' K <- 15  #number of groups
#'
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)
#' X1<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y1 <-X1%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' ## Do cross validation
#' fit.cv <- lsgl::cv(X1, Y1, alpha = 1, lambda = 0.1, intercept = FALSE)
#'
#' ## Cross validation errors (estimated expected generalization error)
#' Err(fit.cv)
#'
#' ## Do the same cross validation using 2 parallel units
#' cl <- makeCluster(2)
#' registerDoParallel(cl)
#'
#' fit.cv <- lsgl::cv(X1, Y1, alpha = 1, lambda = 0.1, intercept = FALSE, use_parallel = TRUE)
#'
#' stopCluster(cl)
#'
#' Err(fit.cv)
#' @author Martin Vincent
#' @importFrom utils packageVersion
#' @importFrom sglOptim sgl_cv
#' @export
cv <- function(x, y,
	intercept = TRUE,
	weights = NULL,
	grouping = NULL,
	groupWeights = NULL,
	parameterWeights = NULL,
	alpha = 1,
	lambda,
	d = 100,
	fold = 10L,
	cv.indices = list(),
	max.threads = NULL,
	use_parallel = FALSE,
	algorithm.config = lsgl.standard.config) {

	# Get call
	cl <- match.call()

	setup <- .process_args(x, y,
		weights = weights,
		intercept = intercept,
		grouping = grouping,
		groupWeights = groupWeights,
		parameterWeights = parameterWeights)

	data <- setup$data

	# Print info
	if(algorithm.config$verbose) {

		n_fold <- if(length(cv.indices) == 0) fold else length(cv.indices)

		cat("\nRunning lsgl", n_fold, "fold cross validation ")
		if(data$sparseX & data$sparseY) {
			cat("(sparse design and response matrices)")
		}
		if(data$sparseX & !data$sparseY) {
			cat("(sparse design matrix)")
		}
		if(!data$sparseX & data$sparseY) {
			cat("(sparse response matrix)")
		}

		cat("\n\n")

		print(data.frame(
			'Samples: ' = print_with_metric_prefix(data$n_samples),
			'Features: ' = print_with_metric_prefix(data$n_covariate),
			'Models: ' = print_with_metric_prefix(ncol(data$data$Y)),
			'Groups: ' = print_with_metric_prefix(length(unique(setup$grouping))),
			'Parameters: ' = print_with_metric_prefix(length(setup$parameterWeights)),
			check.names = FALSE),
		row.names = FALSE, digits = 2, right = TRUE)

		cat("\n")
	}

	res <- sgl_cv(
		module_name = setup$callsym,
		PACKAGE = "lsgl",
		data = data,
		parameterGrouping = setup$grouping,
		groupWeights = setup$groupWeights,
		parameterWeights = setup$parameterWeights,
		alpha =  alpha,
		lambda = lambda,
		d = d,
		fold = fold,
		cv.indices = cv.indices,
		responses = c("link"),
		max.threads = max.threads,
		use_parallel = use_parallel,
		algorithm.config = algorithm.config
		)

	# Add weights
	res$weights <- weights

	res$Yhat <- res$responses$link
	res$responses <- NULL

	res$lsgl_version <- packageVersion("lsgl")
	res$intercept <- intercept
	res$call <- cl

	class(res) <- "lsgl"
	return(res)
}

#' Deprecated cv function
#'
#' @keywords internal
#' @export
lsgl.cv <- function(
  x,
  y,
  intercept = TRUE,
  weights = NULL,
  grouping = NULL,
  groupWeights = NULL,
  parameterWeights = NULL,
  alpha = 1,
  lambda,
  d = 100,
  fold = 10L,
  cv.indices = list(),
  max.threads = NULL,
  use_parallel = FALSE,
  algorithm.config = lsgl.standard.config) {

  warning("lsgl.subsampling is deprecated, use lsgl::subsampling")

  lsgl::cv(
    x = x,
    y = y,
    intercept = intercept,
    weights = weights,
    grouping = grouping,
    groupWeights = groupWeights,
    parameterWeights = parameterWeights,
    alpha = alpha,
    lambda = lambda,
    d = d,
    fold = fold,
    cv.indices = cv.indices,
    max.threads = max.threads,
    use_parallel = use_parallel,
    algorithm.config = algorithm.config
  )

}
