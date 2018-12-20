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

#' @title Compute a lambda sequence for the regularization path
#'
#' @description
#' Computes a decreasing lambda sequence of length \code{d}.
#' The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param y response matrix, matrix of size \eqn{N \times K}.
#' @param intercept should the model include intercept parameters.
#' @param weights sample weights, vector of size \eqn{N \times K}.
#' @param grouping grouping of features, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the feature.
#' @param groupWeights the group weights, a vector of length \eqn{m} (the number of groups).
#' @param parameterWeights a matrix of size \eqn{K \times p}.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param d the length of lambda sequence
#' @param lambda.min the smallest lambda value in the computed sequence.
#' @param lambda.min.rel is lambda.min relative to lambda.max ? (i.e. actual lambda min used is \code{lambda.min*lambda.max}, with \code{lambda.max} the computed maximal lambda value)
#' @param algorithm.config the algorithm configuration to be used.
#' @return a vector of length \code{d} containing the compute lambda sequence.
#' @author Martin Vincent
#' @importFrom methods is
#' @importFrom sglOptim sgl_lambda_sequence
#' @importFrom sglOptim transpose_response_elements
#' @export
lambda <- function(x, y,
	intercept = TRUE,
	weights = NULL,
	grouping = NULL,
	groupWeights = NULL,
	parameterWeights = NULL,
	alpha = 1,
	d = 100L,
	lambda.min,
	lambda.min.rel = FALSE,
	algorithm.config = lsgl.standard.config) {

	setup <- .process_args(x, y,
		weights = weights,
		intercept = intercept,
		grouping = grouping,
		groupWeights = groupWeights,
		parameterWeights = parameterWeights
	)

	data <- setup$data

	lambda <- sgl_lambda_sequence(
		module_name = setup$callsym,
		PACKAGE = "lsgl",
		data = data,
		parameterGrouping = setup$grouping,
		groupWeights = setup$groupWeights,
		parameterWeights = setup$parameterWeights,
		alpha = alpha,
		d = d,
		lambda.min = lambda.min,
		algorithm.config = algorithm.config,
		lambda.min.rel = lambda.min.rel
	)

	return(lambda)
}


#' C interface
#'
#' @keywords internal
#' @export
lsgl_w_xd_yd_sgl_lambda_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  d,
  lambda.min,
  lambda.min.rel,
  algorithm.config) {
  
  .Call(lsgl_w_xd_yd_sgl_lambda, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        d,
        lambda.min,
        lambda.min.rel,
        algorithm.config
  )
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_xd_yd_sgl_lambda_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  d,
  lambda.min,
  lambda.min.rel,
  algorithm.config) {
  
  .Call(lsgl_xd_yd_sgl_lambda, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        d,
        lambda.min,
        lambda.min.rel,
        algorithm.config
  )
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_xs_yd_sgl_lambda_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  d,
  lambda.min,
  lambda.min.rel,
  algorithm.config) {
  
  .Call(lsgl_xs_yd_sgl_lambda, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        d,
        lambda.min,
        lambda.min.rel,
        algorithm.config
  )
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_xd_ys_sgl_lambda_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  d,
  lambda.min,
  lambda.min.rel,
  algorithm.config) {
  
  .Call(lsgl_xd_ys_sgl_lambda, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        d,
        lambda.min,
        lambda.min.rel,
        algorithm.config
  )
}

#' C interface
#'
#' @keywords internal
#' @export
lsgl_xs_ys_sgl_lambda_R <- function(
  data,
  block_dim,
  groupWeights,
  parameterWeights,
  alpha,
  d,
  lambda.min,
  lambda.min.rel,
  algorithm.config) {
  
  .Call(lsgl_xs_ys_sgl_lambda, PACKAGE = "lsgl",
        data,
        block_dim,
        groupWeights,
        parameterWeights,
        alpha,
        d,
        lambda.min,
        lambda.min.rel,
        algorithm.config
  )
}


#' Deprecated lambda function
#'
#' @keywords internal
#' @export
lsgl.lambda <- function(x, y,
  intercept = TRUE,
  weights = NULL,
  grouping = NULL,
  groupWeights = NULL,
  parameterWeights = NULL,
  alpha = 1,
  d = 100L,
  lambda.min,
  lambda.min.rel = FALSE,
  algorithm.config = lsgl.standard.config) {

  warning("lsgl.lambda is deprecated, use lsgl::lambda")

  lsgl::lambda(
    x = x,
    y = y,
    intercept = intercept,
    weights =  weights,
    grouping = grouping,
    groupWeights = groupWeights,
    parameterWeights = parameterWeights,
    alpha = alpha,
    d = d,
    lambda.min  = lambda.min,
    lambda.min.rel = lambda.min.rel,
    algorithm.config = algorithm.config
  )
}
