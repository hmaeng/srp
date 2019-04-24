
#######################################
#   Simulations
#######################################

case1 <-  list(name = "case1",
  Ta = c(0.1,0.2,0.4),
  Tb = srp::truebeta[,1],
  Tm = c(0.01803683)
)

case2 <-  list(name = "case2",
  Ta = c(0.4,-0.5,0.6),
  Tb = srp::truebeta[,2],
  Tm = c(-0.08364838)
)

case3 <-  list(name = "case3",
  Ta = c(0.1,0.2,0.4),
  Tb = srp::truebeta[,3],
  Tm = c(-0.02394883)
)

case4 <-  list(name = "case4",
  Ta = c(0.1,-0.2,0.4),
  Tb = srp::truebeta[,4],
  Tm = c(-0.07421084)
)

sim <- function(N, case, L, maxq, prtn, inisp){

  ### for assesement
  result <- list(qhat=matrix(NA, nrow=N, ncol=8), mspe=matrix(NA, nrow=N, ncol=8), sse=matrix(NA, nrow=N, ncol=8))
  result <- lapply(result, function(x) {colnames(x) <- c("MLR", "FLR", "FLiRTI", "SRPl", "SRPc", "MPDP", "NP", "RDG"); x})

  for(K in 1:N){
    dat <- have.data(case=case, n=300, p=360, K=K, snr=4)

    p <- dim(dat)[1]
    n <- dim(dat)[2]
    n1 <- ceiling(n*prtn)
    n2 <- n-n1

    train_X <- dat[c(1:(p-1)), 1:n1, drop=F]
    test_X <- dat[c(1:(p-1)), (n1+1):n, drop=F]
    train_Y <- dat[p, 1:n1, drop=F]
    test_Y <- dat[p, (n1+1):n, drop=F]

    ##################################################
    ##################### MLR ########################

    mlr.obj <- mlr(dat=dat[, 1:n1, drop=F], maxq=maxq)
    p.mlr <- pred.mlr(mlr.obj, newX=test_X)

    result[[1]][K,1] <- mlr.obj$qhat.mlr
    result[[2]][K,1] <- mean((p.mlr$yhat-c(test_Y))^2, na.rm=T)
    result[[3]][K,1] <- sum((c(case$Tb, case$Ta) - mlr.obj$bhat.mlr)^2)

    ##################################################
    ##################### FLR ########################

    flr.obj <- flr(dat=dat[, 1:n1, drop=F], L=L)
    p.flr <- pred.flr(flr.obj, newX=test_X)

    result[[1]][K,2] <- NA
    result[[2]][K,2] <- mean((p.flr$yhat-c(test_Y))^2, na.rm=T)
    result[[3]][K,2] <- sum((c(case$Tb, case$Ta) - flr.obj$bhat)^2)

    ##################################################
    #################### FLiRTI ######################
    ### available at http://www-bcf.usc.edu/~gareth/research/flrti

    d <- c(2,3,4)
    s <- seq(0.0001, 0.1, length.out=20)
    f.der2 = flrti.cv(Y=c(train_Y), X=t(train_X), sigma=s, deriv=d[1], weight=0.1)
    f.der3 = flrti.cv(Y=c(train_Y), X=t(train_X), sigma=s, deriv=d[2], weight=0.1)
    f.der4 = flrti.cv(Y=c(train_Y), X=t(train_X), sigma=s, deriv=d[3], weight=0.1)

    opt.d = d[which.min(c(f.der2$error[which.min(f.der2$error)],
      f.der3$error[which.min(f.der3$error)],f.der4$error[which.min(f.der4$error)]))]
    min.err=c(which.min(f.der2$error), which.min(f.der3$error), which.min(f.der4$error))
    opt.sigma = s[min.err[opt.d-1]]

    flrti.obj <- flrti(Y=c(train_Y), X=t(train_X), sigma=opt.sigma, deriv=opt.d, weight=0.1, plot=F)
    p.flrti <- predict.flrti(flrti.obj, t(test_X))

    result[[1]][K,3] <- NA
    result[[2]][K,3] <- mean((p.flrti-c(test_Y))^2, na.rm=T)
    result[[3]][K,3] <- sum((c(case$Tb, case$Ta) - c(flrti.obj$beta/(p-1)))^2)

    ##################################################
    ###################### SRPl ######################
    ### use package "srp"

    srpl.obj <- srp.l(dat[, 1:n1, drop=F], maxq=maxq, plot=F)
    p.srpl <-  predict(srpl.obj, x=test_X)

    result[[1]][K,4] <- srpl.obj$qhat
    result[[2]][K,4] <- mean((p.srpl$yhat - c(test_Y))^2, na.rm=T)
    result[[3]][K,4] <- sum((c(case$Tb, case$Ta) - c(srpl.obj$bhat, srpl.obj$ahat))^2)

    ##################################################
    ###################### SRPc ######################
    ### use package "srp"

    srpc.obj <- srp.c(dat[, 1:n1, drop=F], L=L, maxq=maxq, inisp=inisp, plot=F)
    p.srpc <-  predict(srpc.obj, x=test_X)

    result[[1]][K,5] <- srpc.obj$qhat
    result[[2]][K,5] <- mean((p.srpc$yhat - c(test_Y))^2, na.rm=T)
    result[[3]][K,5] <- sum((c(case$Tb, case$Ta) - c(srpc.obj$bhat, srpc.obj$ahat))^2)

    ##################################################
    ##################### MPDP #######################
    ### available at http://www.lsp.ups-tlse.fr/staph/npfda/mpdp-routinesR.txt

    mpdp.obj <- mpdp(t(train_X),c(train_Y), nbmax=5, nbbw=5, Grid=c(1:nrow(train_X)),
      pcvpar=1/3,  kind.of.kernel="quadratic2")
    p.mpdp <- predict.mpdp(mpdp.obj, t(test_X))

    result[[1]][K,6] <- NA
    result[[2]][K,6] <- mean((c(test_Y) - p.mpdp)^2, na.rm=T)
    result[[3]][K,6] <- NA

    ##################################################
    ####################### NP #######################
    ### available at http://www.lsp.ups-tlse.fr/staph/npfda/npfda-routinesR.txt

    np.obj <- funopare.knn.gcv(c(train_Y), t(train_X), t(test_X), 0, nknot=20, c(0,1))
    p.np <- np.obj$Predicted.values

    result[[1]][K,7] <- NA
    result[[2]][K,7] <- mean((c(test_Y) - p.np)^2, na.rm=T)
    result[[3]][K,7] <- NA

    ##################################################
    #################### RDG #########################

    rdg.obj <- rdg(dat=dat[, 1:n1, drop=F])
    p.rdg <- pred.rdg(rdg.obj, newX=test_X)

    result[[1]][K,8] <- NA
    result[[2]][K,8] <- mean((p.rdg$yhat-c(test_Y))^2, na.rm=T)
    result[[3]][K,8] <- sum((c(case$Tb, case$Ta) - rdg.obj$bhat.rdg)^2)

  }

  return(result)

}


sim1 <- sim(N=100, case=case1, L=35, maxq=30, prtn=0.5, inisp=1)
sim2 <- sim(N=100, case=case2, L=35, maxq=30, prtn=0.5, inisp=1)
sim3 <- sim(N=100, case=case3, L=35, maxq=30, prtn=0.5, inisp=1)
sim4 <- sim(N=100, case=case4, L=35, maxq=30, prtn=0.5, inisp=1)


### Table 1
SSE <- rbind(
  paste(round(apply(sim1$sse[,-(6:7),drop=F], 2, mean)*100, 2), "(", round(apply(sim1$sse[,-(6:7),drop=F], 2, sd)*100, 2),")",sep=""),
  paste(round(apply(sim2$sse[,-(6:7),drop=F], 2, mean)*100, 2), "(", round(apply(sim2$sse[,-(6:7),drop=F], 2, sd)*100, 2),")",sep=""),
  paste(round(apply(sim3$sse[,-(6:7),drop=F], 2, mean)*100, 2), "(", round(apply(sim3$sse[,-(6:7),drop=F], 2, sd)*100, 2),")",sep=""),
  paste(round(apply(sim4$sse[,-(6:7),drop=F], 2, mean)*100, 2), "(", round(apply(sim4$sse[,-(6:7),drop=F], 2, sd)*100, 2),")",sep="")
)
colnames(SSE) = c("MLR","FLR","FLiRTI","SRPl", "SRPc", "RDG")
SSE


### Table 2
MSPE <- rbind(
  paste(round(apply(sim1$mspe, 2, mean)*100, 2), "(", round(apply(sim1$mspe, 2, sd)*100, 1),")",sep=""),
  paste(round(apply(sim2$mspe, 2, mean)*100, 2), "(", round(apply(sim2$mspe, 2, sd)*100, 1),")",sep=""),
  paste(round(apply(sim3$mspe, 2, mean)*100, 2), "(", round(apply(sim3$mspe, 2, sd)*100, 1),")",sep=""),
  paste(round(apply(sim4$mspe, 2, mean)*100, 2), "(", round(apply(sim4$mspe, 2, sd)*100, 1),")",sep="")
)
colnames(MSPE) = c("MLR", "FLR", "FLiRTI", "SRPl", "SRPc", "MPDP", "NP", "RDG")
MSPE


### Figure 2
maxfac <- 1+(table(sim1$qhat[,5])==max(table(sim1$qhat[,5])))
barplot(table(sim1$qhat[,5]), col=c("gray","black")[maxfac], cex.names=2, cex.axis=2)
maxfac <- 1+(table(sim2$qhat[,5])==max(table(sim2$qhat[,5])))
barplot(table(sim2$qhat[,5]), col=c("gray","black")[maxfac], cex.names=2, cex.axis=2)
maxfac <- 1+(table(sim3$qhat[,5])==max(table(sim3$qhat[,5])))
barplot(table(sim3$qhat[,5]), col=c("gray","black")[maxfac], cex.names=2, cex.axis=2)
maxfac <- 1+(table(sim4$qhat[,5])==max(table(sim4$qhat[,5])))
barplot(table(sim4$qhat[,5]), col=c("gray","black")[maxfac], cex.names=2, cex.axis=2)


