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

# Check and setup sgl call arguments
#' @keywords internal
#' @importFrom sglOptim add_data
.process_args <- function(x, y,
  weights,
  intercept,
  grouping,
  groupWeights,
  parameterWeights) {

# cast
if( is.null(grouping) )  {
  grouping <- factor(1:ncol(x))
} else {
  grouping <- factor(grouping)
}

# Cat y as matrix
if( is.vector(y) ) {
  y <- matrix(y, nrow = length(y), ncol = 1)
}

# Cat W as matrix
if( is.vector(weights) ) {
  weights <- matrix(weights, nrow = length(weights), ncol = 1)
}

# Validate input
# NOTE ensure ColSums(Y != 0) != 0
if( nrow(x) != nrow(y) ) {
  stop("x and y must have the same number of rows")
}

if(nrow(x) != nrow(y)) {
  stop("x and y must have the same number of rows")
}

if( ! is.null(weights) ) {
  if( ! all(dim(y) == dim(weights))) {
    stop("weights and y must have the same dimensions")
  }
}

# Initialize groupWeights
if( is.null(groupWeights) ) {
  groupWeights <- c(sqrt(ncol(y)*table(grouping)))
}

# Initialize parameterWeights
if( is.null(parameterWeights) ) {
  parameterWeights <-  matrix(1, nrow = ncol(y), ncol = ncol(x))
  dimnames(parameterWeights) <- list(colnames(y), colnames(x))
}

# add intercept
if(intercept) {
  x <- cbind(Intercept = rep(1, nrow(x)), x)
  groupWeights <- c(0, groupWeights)
  parameterWeights <- cbind(rep(0, ncol(y)), parameterWeights)
  grouping <- factor(c("Intercept", as.character(grouping)), levels = c("Intercept", levels(grouping)))
}

# create data
data <- create.sgldata(x, y)
data <- add_data(data, weights, "W")

# Call sglOptim function
callsym <- .get_callsym(data, ! is.null(weights))

setup <- list()
setup$data <- data
setup$callsym <- callsym
setup$grouping <- grouping
setup$groupWeights <- groupWeights
setup$parameterWeights <- parameterWeights

return(setup)

}

# Match with MODULE_NAME in logitsgl.cpp
.get_callsym <- function(data, use_weights) {

  if( use_weights ) {
		obj <- "lsgl_w_"
	} else {
		obj <- "lsgl_"
	}

	return( paste(obj, if(data$sparseX) "xs_" else "xd_", if(data$sparseY) "ys" else "yd", sep = "") )
}
