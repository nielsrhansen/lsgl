#
#     Description of this R script:
#     R tests for linear multiple output sparse group lasso routines.
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


set.seed(100) #  ensures consistency of tests

## Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K

N <- 50 #number of samples
p <- 40 #number of features
K <- 25  #number of groups

B<-matrix(sample(c(rep(1,p*K*0.1),rep(0, p*K-as.integer(p*K*0.1)))),nrow=p,ncol=K)

X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)
Y<-X%*%B+matrix(rnorm(N*K,0,1),N,K)

rownames(X) <- paste("sample", 1:N)
colnames(X) <- paste("feature", 1:p)
rownames(Y) <- paste("sample", 1:N)
colnames(Y) <- paste("model", 1:K)

lambda<-lsgl::lambda(X,Y, alpha=1, lambda.min= 1, intercept=FALSE)

fit <-lsgl::fit(X,Y, alpha=1, lambda = lambda, intercept=FALSE)

## ||B - \beta||_F
if(min(sapply(fit$beta, function(beta) sum((B - beta)^2))) > 11) stop()

# Test predict
res <- predict(fit, X)
stopifnot(all.equal(dimnames(res$Yhat[[10]]), list(rownames(X), colnames(Y))))


print(fit)
features_stat(fit)
parameters_stat(fit)

err <- Err(fit, X)
if(min(err) > 2) stop()
if(min(err) < 0.5) stop()
if(max(err) > 6) stop()

# Test sparse X, Y and intercept

Xsp <- as(X, "CsparseMatrix")
Ysp <- as(Y, "CsparseMatrix")

lambda<-lsgl::lambda(X,Y, alpha=1, lambda.min= 1, intercept = TRUE)
fit <-lsgl::fit(X,Y, alpha=1, lambda = lambda, intercept = TRUE)

lambda<-lsgl::lambda(Xsp,Y, alpha=1, lambda.min= 1, intercept = TRUE)
fit <-lsgl::fit(Xsp,Y, alpha=1, lambda = lambda, intercept = TRUE)

lambda<-lsgl::lambda(X,Ysp, alpha=1, lambda.min= 1, intercept = TRUE)
fit <-lsgl::fit(X,Ysp, alpha=1, lambda = lambda, intercept = TRUE)

lambda<-lsgl::lambda(Xsp,Ysp, alpha=1, lambda.min= 1, intercept = TRUE)
fit <-lsgl::fit(Xsp,Ysp, alpha=1, lambda = lambda, intercept = TRUE)

lambda<-lsgl::lambda(Xsp,Y, alpha=1, lambda.min= 1, intercept=FALSE)
fit <-lsgl::fit(Xsp,Y, alpha=1, lambda = lambda, intercept=FALSE)

lambda<-lsgl::lambda(Xsp,Ysp, alpha=1, lambda.min= 1, intercept=FALSE)
fit <-lsgl::fit(Xsp,Ysp, alpha=1, lambda = lambda, intercept=FALSE)

## Test single fit i.e. K = 1
y <- Y[,1]

lambda<-lsgl::lambda(X,y, alpha=1, lambda.min=1, intercept=FALSE)
fit <-lsgl::fit(X, y, alpha=1, lambda = lambda, intercept=FALSE)
res <- predict(fit, X)

#Check names
stopifnot(all(rownames(res$Yhat) == rownames(X)))

### Navigation tests
print(res)
print(fit)
features_stat(fit)
parameters_stat(fit)

### Test for errors if X or Y contains NA
Xna <- X
Xna[1,1] <- NA

res <- try(lambda<-lsgl::lambda(Xna, Y, alpha=1, lambda.min=.5, intercept=FALSE), silent = TRUE)
if(class(res) != "try-error") stop()

res <- try(fit <-lsgl::fit(Xna, Y, alpha=1, lambda = lambda, intercept=FALSE), silent = TRUE)
if(class(res) != "try-error") stop()

Yna <- Y
Yna[1,1] <- NA

res <- try(lambda<-lsgl::lambda(X, Yna, alpha=1, lambda.min=.5, intercept=FALSE), silent = TRUE)
if(class(res) != "try-error") stop()

res <- try(fit <-lsgl::fit(X, Yna, alpha=1, lambda = lambda, intercept=FALSE), silent = TRUE)
if(class(res) != "try-error") stop()


# test deprecated warnings

assertWarning(
  lambda <- lsgl.lambda(X,Y, alpha=1, lambda.min= 1, intercept = TRUE)
)
assertWarning(
  fit.cv <- lsgl.cv(X,Y, alpha=1, lambda = lambda, intercept = TRUE)
)
