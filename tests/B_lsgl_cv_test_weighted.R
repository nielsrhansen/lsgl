#
#     Description of this R script:
#     R test for linear multiple output using sparse group lasso routines.
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

library(lsgl)
library(methods)

# warnings = errors
options(warn=2)

set.seed(100) # This may be removed, it ensures consistency of the daily tests

## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K

N <- 50 #number of samples
p <- 25 #number of features
K <- 10  #number of groups

B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)
X1<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
Y1 <-X1%*%B+matrix(rnorm(N*K,0,1),N,K)

W <- matrix(1/5, nrow = N, ncol = K)

##Do cross validation
lambda <- lsgl::lambda(X1, Y1, alpha = 1, d = 25, lambda.min = 0.5,  weights = W, intercept = FALSE)
fit.cv <- lsgl::cv(X1, Y1, alpha = 1, lambda = lambda,  weights = W, intercept = FALSE)

## Cross validation errors (estimated expected generalization error)
if(min(Err(fit.cv, loss = "SOVE")) > 0.05) stop()


## Test single fit i.e. K = 1
y <- Y1[,1]
W <- W[,1]

lambda <- lsgl::lambda(X1, y, alpha = 1, d = 25, lambda.min = 0.5,  weights = W, intercept = FALSE)
fit.cv <- lsgl::cv(X1, y, alpha = 1, lambda = lambda,  weights = W, intercept = FALSE)

## Navigation tests
Err(fit.cv)
features_stat(fit.cv)
parameters_stat(fit.cv)
best_model(fit.cv)

### Test for errors if X or Y contains NA
Xna <- X1
Xna[1,1] <- NA

res <- try(fit.cv <- lsgl::cv(Xna, Y1, alpha = 1, lambda = lambda,  weights = W, intercept = FALSE), silent = TRUE)
if(class(res) != "try-error") stop()

Yna <- Y1
Yna[1,1] <- NA

res <- try(fit.cv <- lsgl::cv(X1, Yna, alpha = 1, lambda = lambda,  weights = W, intercept = FALSE), silent = TRUE)
if(class(res) != "try-error") stop()
