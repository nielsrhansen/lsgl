#
#     Description of this R script:
#     R interface for linear multi-response models using sparse group lasso.
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

#' @title Error Rates
#'
#' @description
#' Compute and return an error for each model. The error may be spicifed in the \code{loss} argument.
#'
#' The root-mean-square error (RMSE) is
#' \deqn{\frac{1}{K}\sum_{i = 1}^K \sqrt{\frac{1}{N}\sum_{j=1}^N (Y_{ji}-(X\hat \beta)_{ji})^2}}
#' RMSE is the default error.
#'
#' The objective value error (OVE) is
#' \deqn{\|Y - X\hat \beta\|_F}
#'
#' The scaled objective value error (SOVE) is
#' \deqn{\frac{1}{NK}\|Y - X\hat \beta\|_F}
#'
#' @param object a lsgl object.
#' @param data a design matrix (the \eqn{X} matrix).
#' @param response a matrix of the true responses (the \eqn{Y} matrix).
#' @param loss the loss (error) function. Either a function taking two arguments or
#' one of the following character strings \code{RMSE}, \code{OVE} or \code{SOVE}.
#' @param ... ignored.
#' @return a vector of errors.
#'
#' @author Martin Vincent
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K
#'
#' N <- 100 #number of samples
#' p <- 50 #number of features
#' K <- 15  #number of groups
#'
#' # simulate beta matrix and X matrix
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)
#' X1<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y1 <-X1%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' X2<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y2 <-X2%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' #### Fit models using X1
#' lambda <- lsgl::lambda(X1, Y1, alpha = 1, d = 25L, lambda.min = 5, intercept = FALSE)
#' fit <- lsgl::fit(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE)
#'
#' ## Training errors:
#' Err(fit, X1)
#'
#' ## Errors predicting Y2:
#' Err(fit, X2, Y2)
#'
#' #### Do cross validation
#' fit.cv <- lsgl::cv(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE)
#'
#' ## Cross validation errors (estimated expected generalization error)
#' Err(fit.cv)
#'
#' ## Cross validation errors using objective value error measures
#' Err(fit.cv, loss = "OVE")
#'
#' @importFrom stats predict
#' @importFrom sglOptim Err
#' @importFrom sglOptim compute_error
#' @export
Err.lsgl <- function(object, data = NULL, response = object$Y.true, loss = "RMSE", ... ) {

	if( ! is.function(loss) ) {
		loss <- switch(loss,
			RMSE = function(x,y) mean(sapply(1:length(x), function(i) sqrt(mean((x[[i]] - y[[i]])^2)))),
			OVE = function(x,y) sqrt(sum(sapply(1:length(x), function(i) sum((x[[i]] - y[[i]])^2)))),
			SOVE = function(x,y) 1/(length(x)*length(x[[1]]))*sqrt(sum(sapply(1:length(x), function(i) sum((x[[i]] - y[[i]])^2)))),
			stop("Unknown loss")
		)
	}

	true_response <- response

	if( ! is.null(data) ) {
		object <- predict(object, data)
	}

	return( compute_error(
		x = object,
		response_name = "Yhat",
		true_response = true_response,
		loss = loss)
	)
}

#' @title Nonzero Features
#'
#' @description
#' Extracts the nonzero features for each model.
#'
#' @param object a lsgl object
#' @param ... ignored
#' @return a list of of length \code{nmod(x)} containing the nonzero features (that is nonzero columns of the beta matrices)
#'
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K
#'
#' N <- 100 #number of samples
#' p <- 50 #number of covariates
#' K <- 25  #number of groups
#'
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)
#'
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' lambda<-lsgl::lambda(X,Y, alpha=1, d = 25, lambda.min=.5, intercept=FALSE)
#' fit <-lsgl::fit(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#'
#' # the nonzero features of model 1, 10 and 25
#' features(fit)[c(1,10,25)]
#'
#' # count the number of nonzero features in each model
#' sapply(features(fit), length)
#'
#' @author Martin Vincent
#' @importFrom sglOptim features
#' @export
features.lsgl <- function(object, ...) {

	class(object) <- "sgl" # Use std function

	if( ! is.null(object$beta)) {
		# sgl uses t(beta)
		object$beta <- lapply(object$beta, t)
	}

	return(features(object))
}

#' @title Nonzero Parameters
#'
#' @description
#' Extracts the nonzero parameters for each model.
#'
#' @param object a lsgl object
#' @param ... ignored
#' @return a list of length \code{nmod(x)} containing the nonzero parameters of the models.
#'
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
#'
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' lambda<-lsgl::lambda(X,Y, alpha=1, d = 25, lambda.min=.5, intercept=FALSE)
#' fit <-lsgl::fit(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#'
#' # the nonzero parameters of model 1, 10 and 25
#' parameters(fit)[c(1,10,25)]
#'
#' # count the number of nonzero parameters in each model
#' sapply(parameters(fit), sum)
#'
#' @author Martin Vincent
#' @importFrom sglOptim parameters
#' @export
parameters.lsgl <- function(object, ...) {

	class(object) <- "sgl" # Use std function

	if(!is.null(object$beta)) {
		# sgl uses t(beta)
		object$beta <- lapply(object$beta, t)
	}

	return(parameters(object))
}

#' @title Extract feature statistics
#'
#' @description
#' Extracts the number of nonzero features (or group) in each model.
#'
#' @param object a lsgl object
#' @param ... ignored
#' @return a vector of length \code{nmod(x)} or a matrix containing the number of nonzero features (or group) of the models.
#'
#' @author Martin Vincent
#' @importFrom sglOptim features_stat
#' @export
features_stat.lsgl <- function(object, ...) {

	class(object) <- "sgl" # Use std function

	if(!is.null(object$beta)) {
		# sgl uses t(beta)
		object$beta <- lapply(object$beta, t)
	}

	return(features_stat(object, ...))
}


#' @title Extracting parameter statistics
#'
#' @description
#' Extracts the number of nonzero parameters in each model.
#'
#' @param object a lsgl object
#' @param ... ignored
#' @return a vector of length \code{nmod(x)} or a matrix containing the number of nonzero parameters of the models.
#'
#' @author Martin Vincent
#' @importFrom sglOptim parameters_stat
#' @export
parameters_stat.lsgl <- function(object, ...) {

	class(object) <- "sgl" # Use std function

	if(!is.null(object$beta)) {
		# sgl uses t(beta)
		object$beta <- lapply(object$beta, t)
	}

	return(parameters_stat(object, ...))
}

#' @title Number of Models
#' @description
#' Returns the number of models used for fitting.
#' Note that cv and subsampling objects does not containing any models even though nmod returns a positive number.
#'
#' @param object a lsgl object
#' @param ... ignored
#' @return the number of models in \code{object}
#'
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
#'
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' lambda<-lsgl::lambda(X,Y, alpha=1, d = 25, lambda.min=.5, intercept=FALSE)
#' fit <-lsgl::fit(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#'
#' # the number of models
#' nmod(fit)
#'
#' @author Martin Vincent
#' @importFrom sglOptim nmod
#' @export
nmod.lsgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(nmod(object, ...))
}

#' @title Index of best model
#'
#' @description
#' Returns the index of the best model, in terms of lowest error rate
#' @param object a lsgl object
#' @param ... additional parameters (ignored)
#' @return index of the best model.
#'
#' @author Martin Vincent
#' @importFrom sglOptim best_model
#' @export
best_model.lsgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(best_model(object, "lsgl", ...))
}

#' @title Exstract Fitted Models
#'
#' @description
#' Returns the fitted models, that is the estimated \eqn{\beta} matrices.
#'
#' @param object a lsgl object
#' @param index indices of the models to be returned
#' @param ... ignored
#' @return a list of \eqn{\beta} matrices.
#'
#' @author Martin Vincent
#' @importFrom sglOptim models
#' @export
models.lsgl <- function(object, index = 1:nmod(object), ...) {
	class(object) <- "sgl" # Use std function

	return(models(object, ...))
}

#' @title Nonzero Coefficients
#' @description
#' This function returns the nonzero coefficients (that is the nonzero entries of the \eqn{beta} matrices)
#'
#' @param object a lsgl object
#' @param index indices of the models
#' @param ... ignored
#' @return a list of length \code{length(index)} with nonzero coefficients of the models
#'
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K
#'
#' N <- 100 #number of samples
#' p <- 50 #number of covariates
#' K <- 25  #number of groups
#'
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)
#'
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' lambda<-lsgl::lambda(X,Y, alpha=1, d = 25, lambda.min=.5, intercept=FALSE)
#' fit <-lsgl::fit(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#'
#' # the nonzero coefficients of the models 1, 2 and 5
#' coef(fit, index = c(1,2,5))
#'
#' @author Martin Vincent
#' @importFrom stats coef
#' @export
coef.lsgl <- function(object, index = 1:nmod(object), ...) {

	class(object) <- "sgl" # Use std function

	if(!is.null(object$beta)) {
		# sgl uses t(beta)
		object$beta <- lapply(object$beta, t)
	}

	return(coef(object, index = index, ...))
}


#' Print function for lsgl
#'
#' This function will print some general information about the lsgl object
#'
#' @param x lsgl object
#' @param ... ignored
#'
#' @examples
#'
#' set.seed(100) # This may be removed, it ensures consistency of the daily tests
#'
#' ## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K
#'
#' N <- 100 #number of samples
#' p <- 25 #number of features
#' K <- 15  #number of groups
#'
#' B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)
#'
#' X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
#' Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)
#'
#' lambda<-lsgl::lambda(X,Y, alpha=1, d = 25, lambda.min= 5, intercept=FALSE)
#' fit <-lsgl::fit(X,Y, alpha=1, lambda = lambda, intercept=FALSE)
#'
#' # Print some information about the estimated models
#' fit
#'
#' ## Cross validation
#' fit.cv <- lsgl::cv(X, Y, alpha = 1, lambda = lambda, intercept = FALSE)
#'
#'# Print some information
#' fit.cv
#'
#' @author Martin Vincent
#' @importFrom sglOptim sgl_print
#' @method print lsgl
#' @export
print.lsgl <- function(x, ...) {
	sgl_print(x)
}
