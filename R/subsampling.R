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

#' @title Subsampling
#' @description
#' Linear multiple output subsampling using multiple possessors
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param y response matrix, matrix of size \eqn{N \times K}.
#' @param intercept should the model include intercept parameters.
#' @param weights sample weights, vector of size \eqn{N \times K}.
#' @param grouping grouping of features, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the feature.
#' @param groupWeights the group weights, a vector of length \eqn{m} (the number of groups).
#' @param parameterWeights a matrix of size \eqn{K \times p}.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda lambda.min relative to lambda.max or the lambda sequence for the regularization path (that is a vector or a list of vectors with the lambda sequence for the subsamples).
#' @param d length of lambda sequence (ignored if \code{length(lambda) > 1})
#' @param train a list of training samples, each item of the list corresponding to a subsample.
#' Each item in the list must be a vector with the indices of the training samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{test} list.
#' @param test a list of test samples, each item of the list corresponding to a subsample.
#' Each item in the list must be vector with the indices of the test samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{training} list.
#' @param collapse if \code{TRUE} the results for each subsample will be collapse into one result (this is useful if the subsamples are not overlapping)
#' @param max.threads Deprecated (will be removed in 2018),
#' instead use \code{use_parallel = TRUE} and registre parallel backend (see package 'doParallel').
#' The maximal number of threads to be used.
#' @param use_parallel If \code{TRUE} the \code{foreach} loop will use \code{\%dopar\%}. The user must registre the parallel backend.
#' @param algorithm.config the algorithm configuration to be used.
#' @return
#' \item{Yhat}{if \code{collapse = FALSE} then a list of length \code{length(test)} containing the predicted responses for each of the test sets. If \code{collapse = TRUE} a list of length \code{length(lambda)}}
#' \item{Y.true}{a list of length \code{length(test)} containing the true responses of the test samples}
#' \item{features}{number of features used in the models}
#' \item{parameters}{number of parameters used in the models.}
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K
#'
#' N <- 100 #number of samples
#' p <- 50 #number of features
#' K <- 25  #number of groups
#'
#' B <- matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)
#' X1 <- matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y1 <- X1%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' ## Do cross subsampling
#'
#' train <- replicate(2, sample(1:N, 50), simplify = FALSE)
#' test <- lapply(train, function(idx) (1:N)[-idx])
#'
#' lambda <- lapply(train, function(idx)
#'	lsgl::lambda(
#'		x = X1[idx,],
#'		y = Y1[idx,],
#'		alpha = 1,
#'		d = 15L,
#'		lambda.min = 5,
#'		intercept = FALSE)
#'	)
#'
#' fit.sub <- lsgl::subsampling(
#'	 x = X1,
#'	 y = Y1,
#'	 alpha = 1,
#'	 lambda = lambda,
#'	 train = train,
#'	 test = test,
#'	 intercept = FALSE
#' )
#'
#' Err(fit.sub)
#'
#' ## Do the same cross subsampling using 2 parallel units
#' cl <- makeCluster(2)
#' registerDoParallel(cl)
#'
#' # Run subsampling
#' # Using a lambda sequence ranging from the maximal lambda to 0.1 * maximal lambda
#' fit.sub <- lsgl::subsampling(
#'	 x = X1,
#'	 y = Y1,
#'	 alpha = 1,
#'	 lambda = 0.1,
#'	 train = train,
#'	 test = test,
#'	 intercept = FALSE
#' )
#'
#' stopCluster(cl)
#'
#' Err(fit.sub)
#' @author Martin Vincent
#' @importFrom utils packageVersion
#' @importFrom sglOptim sgl_subsampling
#' @importFrom methods is
#' @export
subsampling <- function(
  x,
  y,
  intercept = TRUE,
  weights = NULL,
  grouping = NULL,
  groupWeights = NULL,
  parameterWeights =  NULL,
  alpha = 1,
  lambda,
  d = 100,
  train,
  test,
  collapse = FALSE,
  max.threads = NULL,
  use_parallel = FALSE,
  algorithm.config = lsgl.standard.config)
{

	# Get call
	cl <- match.call()

	setup <- .process_args(x, y,
		weights = weights,
		intercept = intercept,
		grouping = grouping,
		groupWeights = groupWeights,
		parameterWeights = parameterWeights
	)

	data <- setup$data

	# Print info
	if(algorithm.config$verbose) {

		cat("\nRunning lsgl subsampling with ", length(train)," subsamples ")
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

	res <- sgl_subsampling(
		module_name = setup$callsym,
		PACKAGE = "lsgl",
		data = data,
		parameterGrouping = setup$grouping,
		groupWeights = setup$groupWeights,
		parameterWeights = setup$parameterWeights,
		alpha =  alpha,
		lambda = lambda,
		d = d,
		training = train,
		test = test,
		collapse = collapse,
		responses = "link",
		max.threads = max.threads,
		use_parallel = use_parallel,
		algorithm.config = algorithm.config
	)

	# Add weights
	res$weights <- weights

	# Responses
	res$Yhat <- res$responses$link
	res$responses <- NULL

	res$lsgl_version <- packageVersion("lsgl")
	res$intercept <- intercept
	res$call <- cl

	class(res) <- "lsgl"
	return(res)
}

#' Deprecated subsampling function
#'
#' @keywords internal
#' @export
lsgl.subsampling <- function(
  x,
  y,
  intercept = TRUE,
  weights = NULL,
  grouping = NULL,
  groupWeights = NULL,
  parameterWeights =  NULL,
  alpha = 1,
  lambda,
  d = 100,
  train,
  test,
  collapse = FALSE,
  max.threads = NULL,
  use_parallel = FALSE,
  algorithm.config = lsgl.standard.config) {

  warning("lsgl.subsampling is deprecated, use lsgl::subsampling")

  lsgl::subsampling(
    x = x,
    y = y,
    intercept = intercept,
    weights = weights,
    grouping = grouping,
    groupWeights = groupWeights,
    parameterWeights =  parameterWeights,
    alpha = alpha,
    lambda = lambda,
    d = d,
    train = train,
    test = test,
    collapse = collapse,
    max.threads = max.threads,
    use_parallel = use_parallel,
    algorithm.config = algorithm.config
  )
}
