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
library(tools)
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

rownames(X1) <- paste("sample", 1:N)
colnames(X1) <- paste("feature", 1:p)
rownames(Y1) <- paste("sample", 1:N)
colnames(Y1) <- paste("model", 1:K)

##Do cross validation
lambda <- lsgl::lambda(X1, Y1, alpha = 1, d = 25, lambda.min = 0.5, intercept = FALSE)

train <- replicate(2, 1:N, simplify = FALSE)
test <- list(1:20, 21:N)

fit.sub <- lsgl::subsampling(X1, Y1, alpha = 1, lambda = lambda, intercept = FALSE, train = train, test = test)

stopifnot(all.equal(dimnames(fit.sub$Yhat[[1]][[10]]), list(rownames(X1)[test[[1]]], colnames(Y1))))
stopifnot(all.equal(dimnames(fit.sub$Yhat[[2]][[10]]), list(rownames(X1)[test[[2]]], colnames(Y1))))

## Navigation tests
Err(fit.sub)
features_stat(fit.sub)
parameters_stat(fit.sub)
best_model(fit.sub)

Xsp <- as(X1, "CsparseMatrix")
Ysp <- as(Y1, "CsparseMatrix")

fit.sub <- lsgl::subsampling(X1, Y1, alpha = 1, lambda = lambda, intercept = TRUE, train = train, test = test)
fit.sub <- lsgl::subsampling(Xsp, Y1, alpha = 1, lambda = lambda, intercept = TRUE, train = train, test = test)
fit.sub <- lsgl::subsampling(X1, Ysp, alpha = 1, lambda = lambda, intercept = TRUE, train = train, test = test)
fit.sub <- lsgl::subsampling(Xsp, Ysp, alpha = 1, lambda = lambda, intercept = FALSE, train = train, test = test)
fit.sub <- lsgl::subsampling(Xsp, Y1, alpha = 1, lambda = lambda, intercept = FALSE, train = train, test = test)
fit.sub <- lsgl::subsampling(Xsp, Ysp, alpha = 1, lambda = lambda, intercept = FALSE, train = train, test = test)


## Test single fit i.e. K = 1
y <- Y1[,1]
fit.sub <- lsgl::subsampling(X1, y, alpha = 1, lambda = lambda, intercept = TRUE, train = train, test = test)


# test deprecated warnings

assertWarning(
  fit.cv <- lsgl.subsampling(X1, y, alpha = 1, lambda = lambda, intercept = TRUE, train = train, test = test)
)
