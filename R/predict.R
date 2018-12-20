#
#     Description of this R script:
#     R interface for linear multi-response sparse group lasso routines.
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

#' @title Predict
#'
#' @description
#' Compute the predicted response matrix for a new data set.
#'
#' @param object an object of class lsgl, produced with \code{lsgl}.
#' @param x a data matrix of size \eqn{N_\textrm{new} \times p}.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param ... ignored.
#' @return
#' \item{Yhat}{the predicted response matrix (of size \eqn{N_\textrm{new} \times K})}
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
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)
#' X1<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y1 <-X1%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' X2<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y2 <-X2%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' #### Fit models using X1
#' lambda <- lsgl::lambda(X1, Y1, alpha = 1, d = 25L, lambda.min = 0.5, intercept = FALSE)
#' fit <- lsgl::fit(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE)
#'
#' # Predict Y2 using the estimated models and X2
#' res <- predict(fit, X2)
#'
#'  # Frobenius norm \|Y2hat - Y2\|_F
#' sapply(res$Yhat, function(y) sqrt(sum((y - Y2)^2)))
#'
#' # This is proportional to the errors compute with Err:
#' Err(fit, X2, Y2)*length(Y2)
#'
#' @author Martin Vincent
#' @method predict lsgl
#' @importFrom methods is
#' @importFrom methods as
#' @importFrom sglOptim sgl_predict
#' @importFrom sglOptim create.sgldata
#' @export
predict.lsgl <- function(object, x, sparse.data = is(x, "sparseMatrix"), ...)
{

	# Get call
	cl <- match.call()

	if(is.null(object$beta)) stop("No models found -- missing beta")

	object$beta <- lapply(object$beta, t)

	if(object$intercept){
		# add intercept
		x <- cBind(Intercept = rep(1, nrow(x)), x)
	}

	#Check dimension of x
	if(dim(object$beta[[1]])[2] != ncol(x)) stop("x has wrong dimension")

	data <- create.sgldata(
		x = x,
		y = NULL,
		response_dimension = nrow(object$beta[[1]]),
		response_names = rownames(object$beta[[1]]),
		sparseX = sparse.data,
		sparseY = FALSE
	)

	res <- sgl_predict(
		module_name = "lsgl_xd_yd",
		PACKAGE = "lsgl",
		object = object,
		data = data,
		responses = "link")

	#Responses
	res$Yhat <-res$responses$link
	res$responses <- NULL

	res$call <- cl
  res$lsgl_version <- packageVersion("lsgl")
	class(res) <- "lsgl"

	return(res)
}


#' C interface
#'
#' @keywords internal
#' @export
lsgl_xd_yd_sgl_predict_R <- function(data, beta) {
  .Call(lsgl_xd_yd_sgl_predict, PACKAGE = "lsgl", data, beta)
}


#' C interface
#'
#' @keywords internal
#' @export
lsgl_xs_yd_sgl_predict_R <- function(data, beta) {
  .Call(lsgl_xs_yd_sgl_predict, PACKAGE = "lsgl", data, beta)
}


#' C interface
#'
#' @keywords internal
#' @export
lsgl_xd_ys_sgl_predict_R <- function(data, beta) {
  .Call(lsgl_xd_ys_sgl_predict, PACKAGE = "lsgl", data, beta)
}


#' C interface
#'
#' @keywords internal
#' @export
lsgl_xs_ys_sgl_predict_R <- function(data, beta) {
  .Call(lsgl_xs_ys_sgl_predict, PACKAGE = "lsgl", data, beta)
}

