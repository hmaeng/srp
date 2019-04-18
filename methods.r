library(srp)
library(fda)
library(mgcv)
library(lpSolve)
library(tidyr)
library(MASS)
library(xts)
source("flrti.txt") # http://www-bcf.usc.edu/~gareth/research/flrti
file("http://www.lsp.ups-tlse.fr/staph/npfda/npfda-routinesR.txt", open="r")
source("http://www.lsp.ups-tlse.fr/staph/npfda/npfda-routinesR.txt")
file("http://www.lsp.ups-tlse.fr/staph/npfda/mpdp-routinesR.txt",open="r")
source("http://www.lsp.ups-tlse.fr/staph/npfda/mpdp-routinesR.txt")


have.data <- function(case, n, p, K, snr){
  dat0 = matrix(NA, ncol=n, nrow=p)
  for(i in 1:n){
    set.seed(n*(K-1)+i)
    dat0[,i] = arima.sim(model=list(ar=rev(case$Ta)), n=p)
  }

  q0 = length(case$Ta)
  m = p-q0-1
  XV = t(dat0)[,-p]
  XB = XV[ ,1:m,drop=F]%*%matrix(case$Tb, ncol=1) + XV[ ,(m+1):(m+q0),drop=F]%*%matrix(case$Ta, ncol=1)
  errvar = var(XB)/snr
  set.seed(30000+K)
  err = matrix(rnorm(n, mean=0, sd=sqrt(errvar)), ncol=1)

  Y = matrix(matrix(case$Tm, nrow=n, ncol=1)+XB+err, nrow=n, ncol=1)
  dat = t(cbind(XV, Y))

  return(dat)
}

flr <- function(dat, L=35){

  p <- dim(dat)[1]
  n <- dim(dat)[2]

  X <- dat[c(1:(p-1)),,drop=F]
  Y <- dat[p,,drop=F]
  cY <- Y-mean(Y)

  t <- seq(0, 1, length.out=p-1)*(p-1)

  B.basis <- create.bspline.basis(rangeval=c(0, p-1), norder=4, nbasis=L)

  Xfd = smooth.basis(t, X, B.basis)$fd  ## functional data object for train_X
  cXfd = center.fd(Xfd)

  J <- inprod(B.basis, B.basis)
  R <- eval.penalty(B.basis, int2Lfd(2))
  W <- t(cXfd$coefs)%*%J

  mgfit.flr <- magic(y=c(cY), X=W, sp=1, S=list(R), off=1, gcv=T)
  etahat.flr <- mgfit.flr$b
  bhat.flr <- eval.basis(t, B.basis)%*%etahat.flr
  muhat.flr <- mean(Y)-c(t(rowMeans(mean(Xfd)$coefs))%*%J%*%matrix(etahat.flr, ncol=1))

  results <- list("muhat.flr" = muhat.flr, "bhat.flr" = bhat.flr, "etahat.flr" = etahat.flr, "J" = J, "L" = L)

  class(results) <- "flr"
  return(results)
}

pred.flr <- function(object, newX){

  p <- dim(newX)[1] + 1
  n <- dim(newX)[2]
  t <- seq(0, 1, length.out=p-1)*(p-1)

  B.basis <- create.bspline.basis(rangeval=c(0, p-1), norder=4, nbasis=object$L)
  Xfd <- smooth.basis(t, newX, B.basis)$fd
  W <- t(Xfd$coefs)%*%(object$J)
  yhat <- W%*%object$etahat.flr + object$muhat.flr

  return(list("yhat" = yhat))
}

mlr <- function(dat, maxq){
  p <- dim(dat)[1]
  n <- dim(dat)[2]

  Y <- dat[p,]
  cY <- Y-mean(Y)

  SICq.mlr <- rep(0, maxq)

  for(q in 1:maxq){

    X <- dat[(p-q):(p-1),,drop=F]
    cX <- X-apply(X, 1, mean)

    linfit <- lm(cY ~ t(cX)-1, data=data.frame(cbind(matrix(cY, ncol=1), t(cX))))
    muhat.mlr <- mean(Y) - c(matrix(coef(linfit), nrow=1)%*%matrix(apply(X, 1, mean), ncol=1))

    SICq.mlr[q] <- log(mean((Y-(t(X)%*%matrix(coef(linfit), ncol=1) + muhat.mlr))^2))*n + log(n)*(q+1)
  }
  qhat.mlr <- which.min(SICq.mlr)

  X <- dat[(p-qhat.mlr):(p-1),,drop=F]
  cX <- X-apply(X, 1, mean)
  linfit <- lm(cY ~ t(cX)-1, data=data.frame(cbind(cY, t(cX))))
  bhat.mlr <- c(rep(0, p-qhat.mlr-1), coef(linfit))
  muhat.mlr <- mean(Y) - c(matrix(coef(linfit), nrow=1)%*%matrix(apply(X, 1, mean), ncol=1))

  results <- list("muhat.mlr" = muhat.mlr, "coeflinfit"=coef(linfit), "bhat.mlr" = bhat.mlr, "qhat.mlr" = qhat.mlr)

  class(results) <- "mlr"
  return(results)
}

pred.mlr <- function(object, newX){
  p <- dim(newX)[1] + 1
  X <- newX[(p-object$qhat.mlr):(p-1),,drop=F]
  yhat <- t(X)%*%matrix(object$coeflinfit, ncol=1) + object$muhat.mlr
  return(list("yhat" = yhat))
}

rdg <- function(dat){
  p <- dim(dat)[1]
  n <- dim(dat)[2]

  X <- dat[1:(p-1),,drop=F]
  cX <- X-apply(X, 1, mean)
  Y <- dat[p,]
  cY <- Y-mean(Y)

  mgfit.rdg <- magic(y=c(cY), X=t(cX), sp=c(1), S=list(diag(p-1)), off=c(1), gcv=T)
  bhat.rdg <-  mgfit.rdg$b
  muhat.rdg <- mean(Y) - c(matrix(bhat.rdg, nrow=1)%*%matrix(apply(X, 1, mean), ncol=1))

  results <- list("muhat.rdg" = muhat.rdg, "bhat.rdg" = bhat.rdg)

  class(results) <- "rdg"
  return(results)
}

pred.rdg <- function(object, newX){
  yhat <- t(newX)%*%matrix(object$bhat.rdg, ncol=1) + object$muhat.rdg
  return(list("yhat" = yhat))
}

### downloaded from http://www.lsp.ups-tlse.fr/staph/npfda/mpdp-routinesR.txt
### and slightly adjusted as the Bwdselection does not work when only one variable is added in Fwdselection
predict.mpdp <- function(object, CURVPRED)
{
  if(is.vector(CURVPRED)) CURVPRED <- as.matrix(CURVPRED)
  ### Adjustment: the Bwdselection does not work when only one variable is added in Fwdselection
  if(is.numeric(object$Bwdselection)==F){
    VECEST <- object$predictors[,object$Fwdselection, drop=F]
    VECPRED <- CURVPRED[,object$Fwdselection, drop=F]
  }else{
    VECEST <- object$predictors[,object$Bwdselection]
    VECPRED <- CURVPRED[,object$Bwdselection]
  }
  Meancurve <- apply(VECEST,2,mean)
  CCOVARIATES <- t(t(VECEST)-Meancurve)
  Stdev <- sqrt(apply(CCOVARIATES^2,2,mean))
  COVARIATE <- t(t(CCOVARIATES)/Stdev)
  CVECPRED <- t(t(VECPRED)-Meancurve)
  PRED <- t(t(CVECPRED)/Stdev)
  nind <- nrow(COVARIATE)
  kernel <- get(object$kind.of.kernel)
  nind2 <- nrow(PRED)
  Ind1 <- rep(1:nind2,nind)
  Ind2 <- rep(1:nind,rep(nind2,nind))
  COVAR.REP.ROW <- COVARIATE[Ind2,,drop=F]
  PRED.ALT.ROW <- PRED[Ind1,,drop=F]
  D2 <- (PRED.ALT.ROW-COVAR.REP.ROW)^2
  D2OVERBW <- D2/(object$bandwidth^2)
  Dist <- sqrt(apply(D2OVERBW,1,sum))
  Ker <- rep(0,nind2*nind)
  Ker[Dist>0] <- kernel(Dist[Dist>0])
  KERNEL <- matrix(Ker, nrow=nind2)
  Denom <- apply(KERNEL,1,sum)
  NUM <- KERNEL %*% COVARIATE
  XBAR <- NUM/Denom
  Ynum <- as.vector(KERNEL %*% object$responses)
  Ybar <- Ynum/Denom
  TP1 <- crossprod(t(KERNEL)*object$responses,COVARIATE)
  TP2 <- XBAR*Ynum
  TP1M2 <- TP1-TP2
  Predictions <- 0
  for(ii in 1:nind2){
    CCOVAR <- t(t(COVARIATE)-XBAR[ii,])
    TP3 <- crossprod(CCOVAR,(CCOVAR*KERNEL[ii,]))
    ### Adjustment: NA is returned due to a numerical issue
    if(is.na(abs(det(TP3))) || abs(det(TP3))<1e-8){
      Predictions[ii] <- Ybar[ii]
    }else{
      Bb <- solve(TP3,TP1M2[ii,])
      Cpred <- t(t(PRED)-XBAR[ii,])
      Predictions[ii] <- Ybar[ii]+crossprod(Bb,Cpred[ii,])
    }
  }
  return(Predictions)
}

### downloaded from http://www.lsp.ups-tlse.fr/staph/npfda/npfda-routinesR.txt
### and slightly adjusted as error occurs when the sample size is smaller than 20
funopare.knn.gcv <- function(Response, CURVES, PRED, ..., kind.of.kernel = "quadratic", semimetric = "deriv"){
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr")
    SEMIMETRIC1 <- sm(Response, CURVES, CURVES, ...)
  else SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  kernel <- get(kind.of.kernel)
  n1 <- ncol(SEMIMETRIC1)
  step <- ceiling(n1/100)
  if(step == 0)
    step <- 1
  ### ERR occurs in "Knearest <- seq(from = 10, to = n1 %/% 2, by = step)" when n1/2 is smaller than 10
  #Knearest <- seq(from = 10, to = n1 %/% 2, by = step)
  Knearest <- seq(from = 1, to = n1 %/% 2, by = step) # adjusted
  kmax <- max(Knearest)
  # the vector Knearest contains the sequence of the
  # k-nearest neighbours used for computing the optimal bandwidth
  Response.estimated <- 0
  Bandwidth.opt <- 0
  HAT.RESP <- matrix(0, nrow = n1, ncol = length(Knearest))
  BANDWIDTH <- matrix(0, nrow = n1, ncol = kmax)
  for(i in 1:n1) {
    Norm.diff <- SEMIMETRIC1[, i]
    # "norm.order" gives the sequence k_1, k_2,... such that
    # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
    Norm.order <- order(Norm.diff)
    # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),...,
    # dq(X_{j_{kamx+2}},X_i)
    zz <- sort(Norm.diff)[2:(kmax + 2)]
    # BANDWIDTH[i, l-1] contains (dq(X_{j_l},X_i) +
    # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
    BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
    z <- zz[ - (kmax + 1)]
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    UMAT <- ZMAT/BANDWIDTH[i,  ]
    KMAT <- kernel(UMAT)
    KMAT[col(KMAT) > row(KMAT)] <- 0
    Ind.curves <- Norm.order[2:(kmax + 1)]
    Ind.resp <- Response[Ind.curves]
    YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
    HAT.RESP[i,  ] <- apply(YMAT[Knearest,  ] * KMAT[Knearest,  ],
      1, sum)/apply(KMAT[Knearest,  ], 1, sum)
  }
  CRITERIUM <- (HAT.RESP - Response)^2
  Criterium <- apply(CRITERIUM, 2, sum)
  index.opt <- order(Criterium)[1]
  Response.estimated <- HAT.RESP[, index.opt]
  knearest.opt <- Knearest[index.opt]
  Bandwidth.opt <- BANDWIDTH[, knearest.opt]
  Mse.estimated <- sum((Response.estimated - Response)^2)/n1
  if(twodatasets) {
    if(semimetric == "mplsr")
      SEMIMETRIC2 <- sm(Response, CURVES, PRED, ...)
    else SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    Bandwidth2 <- 0
    n2 <- ncol(SEMIMETRIC2)
    for(k in 1:n2) {
      Sm2k <- SEMIMETRIC2[, k]
      Bandwidth2[k] <- sum(sort(Sm2k)[knearest.opt:(knearest.opt+1)])*0.5
    }
    KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    Denom <- apply(KERNEL, 2, sum)
    RESPKERNEL <- KERNEL * Response
    Response.predicted <- apply(RESPKERNEL, 2, sum)/Denom
    return(list(Estimated.values = Response.estimated,
      Predicted.values = Response.predicted, Bandwidths =
        Bandwidth.opt, knearest.opt = knearest.opt, Mse =
        Mse.estimated))
  }else {
    return(list(Estimated.values = Response.estimated, Bandwidths
      = Bandwidth.opt, knearest.opt = knearest.opt, Mse =
        Mse.estimated))
  }
}

