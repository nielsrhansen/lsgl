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

#' @title Fit a linear multiple output model using sparse group lasso
#'
#' @description
#' For a linear multiple output model with \eqn{p} features (covariates) dived into \eqn{m} groups using sparse group lasso.
#'
#' @details
#' This function computes a sequence of minimizers (one for each lambda given in the \code{lambda} argument) of
#' \deqn{\frac{1}{N}\|Y-X\beta\|_F^2 + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
#' where \eqn{\|\cdot\|_F} is the frobenius norm.
#' The vector \eqn{\beta^{(J)}} denotes the parameters associated with the \eqn{J}'th group of features.
#' The group weights are denoted by \eqn{\gamma \in [0,\infty)^m} and the parameter weights by \eqn{\xi \in [0,\infty)^n}.
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param y response matrix, matrix of size \eqn{N \times K}.
#' @param intercept should the model include intercept parameters.
#' @param weights sample weights, vector of size \eqn{N \times K}.
#' @param grouping grouping of features, a factor or vector of length \eqn{p}.
#' Each element of the factor/vector specifying the group of the feature.
#' @param groupWeights the group weights, a vector of length \eqn{m} (the number of groups).
#' @param parameterWeights a matrix of size \eqn{K \times p}.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda lambda.min relative to lambda.max or the lambda sequence for the regularization path.
#' @param d length of lambda sequence (ignored if \code{length(lambda) > 1})
#' @param algorithm.config the algorithm configuration to be used.
#'
#' @return
#' \item{beta}{the fitted parameters -- the list \eqn{\hat\beta(\lambda(1)), \dots, \hat\beta(\lambda(d))} of length \code{length(return)}.
#' With each entry of list holding the fitted parameters, that is matrices of size \eqn{K\times p} (if \code{intercept = TRUE} matrices of size \eqn{K\times (p+1)})}
#' \item{loss}{the values of the loss function.}
#' \item{objective}{the values of the objective function (i.e. loss + penalty).}
#' \item{lambda}{the lambda values used.}
#' @author Martin Vincent
#' @examples
#'
#' set.seed(100) # This may be removed, ensures consistency of tests
#'
#' # Simulate from Y = XB + E,
#' # the dimension of Y is N x K, X is N x p, B is p x K
#'
#' N <- 50 # number of samples
#' p <- 50 # number of features
#' K <- 25 # number of groups
#'
#' B <- matrix(
#'	sample(c(rep(1,p*K*0.1), rep(0, p*K-as.integer(p*K*0.1)))),
#' 	nrow = p,ncol = K)
#'
#' X <- matrix(rnorm(N*p,1,1), nrow=N, ncol=p)
#' Y <- X %*% B + matrix(rnorm(N*K,0,1), N, K)
#'
#' fit <-lsgl::fit(X,Y, alpha=1, lambda = 0.1, intercept=FALSE)
#'
#' ## ||B - \beta||_F
#' sapply(fit$beta, function(beta) sum((B - beta)^2))
#'
#' ## Plot
#' par(mfrow = c(3,1))
#' image(B, main = "True B")
#' image(
#'	x = as.matrix(fit$beta[[100]]),
#' 	main = paste("Lasso estimate (lambda =", round(fit$lambda[100], 2), ")")
#' )
#' image(solve(t(X)%*%X)%*%t(X)%*%Y, main = "Least squares estimate")
#'
#' # The training error of the models
#' Err(fit, X, loss="OVE")
#' # This is simply the loss function
#' sqrt(N*fit$loss)
#'
#' @import Matrix
#' @importFrom utils packageVersion
#' @importFrom methods is
#' @importFrom sglOptim sgl_fit
#' @importFrom sglOptim print_with_metric_prefix
#' @export
fit <- function(x, y,
		intercept = TRUE,
		weights = NULL,
		grouping = NULL,
		groupWeights = NULL,
		parameterWeights = NULL,
		alpha = 1,
		lambda,
		d = 100,
		algorithm.config = lsgl.standard.config)
{
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

		cat("\nRunning lsgl ")
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
			row.names = FALSE, digits = 2, right = TRUE
		)

		cat("\n")
	}


	# Call sglOptim
	res <- sgl_fit(
		module_name = setup$callsym,
		PACKAGE = "lsgl",
		data = data,
		parameterGrouping = setup$grouping,
		groupWeights = setup$groupWeights,
		parameterWeights = setup$parameterWeights,
		alpha = alpha,
		lambda = lambda,
		d = d,
		algorithm.config = algorithm.config
	)

	# Add weights
	res$weights <- weights

# Transpose all beta matrices
	res$beta <- lapply(res$beta, t)

	res$intercept <- intercept
	res$lsgl_version <- packageVersion("lsgl")
	res$call <- cl

	class(res) <- "lsgl"
	return(res)
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_w_xd_yd_sgl_fit_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  lambda,
  idx,
  algorithm.config) {
  
  .Call(lsgl_w_xd_yd_sgl_fit, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        lambda,
        idx,
        algorithm.config
  )
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_xd_yd_sgl_fit_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  lambda,
  idx,
  algorithm.config) {
  
  .Call(lsgl_xd_yd_sgl_fit, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        lambda,
        idx,
        algorithm.config
  )
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_xs_yd_sgl_fit_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  lambda,
  idx,
  algorithm.config) {
  
  .Call(lsgl_xs_yd_sgl_fit, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        lambda,
        idx,
        algorithm.config
  )
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_xd_ys_sgl_fit_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  lambda,
  idx,
  algorithm.config) {
  
  .Call(lsgl_xd_ys_sgl_fit, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        lambda,
        idx,
        algorithm.config
  )
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_xs_ys_sgl_fit_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  lambda,
  idx,
  algorithm.config) {
  
  .Call(lsgl_xs_ys_sgl_fit, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        lambda,
        idx,
        algorithm.config
  )
}



#' Deprecated fit function
#'
#' @keywords internal
#' @export
lsgl <- function(
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
  algorithm.config = lsgl.standard.config) {

  warning("lsgl is deprecated, use lsgl::fit")

  lsgl::fit(
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
    algorithm.config = algorithm.config
  )
}

