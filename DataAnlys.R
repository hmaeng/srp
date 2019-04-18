#######################################
#   Data applications
#######################################

DATAanlys <- function(dat0, N, L, maxq, prtn, nknot, nbmax, nbbw, pcvpar, s, inisp){

  ### for assesement
  result <- list(qhat=matrix(NA, nrow=N, ncol=8), mspe=matrix(NA, nrow=N, ncol=8))
  result <- lapply(result, function(x) {colnames(x) <- c("MLR", "FLR", "FLiRTI", "SRPl", "SRPc", "MPDP", "NP", "RDG"); x})

  for(K in 1:N){

    p <- dim(dat0)[1]
    n <- dim(dat0)[2]

    set.seed(K)
    samp <- sample(n)
    dat <- dat0[,c(samp)]

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

    ##################################################
    ##################### FLR ########################

    flr.obj <- flr(dat=dat[, 1:n1, drop=F], L=L)
    p.flr <- pred.flr(flr.obj, newX=test_X)

    result[[1]][K,2] <- NA
    result[[2]][K,2] <- mean((p.flr$yhat-c(test_Y))^2, na.rm=T)

    ##################################################
    #################### FLiRTI ######################
    ### available at http://www-bcf.usc.edu/~gareth/research/flrti

    d <- c(2,3,4)
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

    ##################################################
    ###################### SRPl ######################
    ### use package "srp"

    srpl.obj <- srp.l(dat[, 1:n1, drop=F], maxq=maxq, plot=F)
    p.srpl <-  predict(srpl.obj, x=test_X)

    result[[1]][K,4] <- srpl.obj$qhat
    result[[2]][K,4] <- mean((p.srpl$yhat - c(test_Y))^2, na.rm=T)

    ##################################################
    ###################### SRPc ######################
    ### use package "srp"

    srpc.obj <- srp.c(dat[, 1:n1, drop=F], L=L, maxq=maxq, inisp=inisp, plot=F)
    p.srpc <-  predict(srpc.obj, x=test_X)

    result[[1]][K,5] <- srpc.obj$qhat
    result[[2]][K,5] <- mean((p.srpc$yhat - c(test_Y))^2, na.rm=T)

    ##################################################
    ##################### MPDP #######################
    ### available at http://www.lsp.ups-tlse.fr/staph/npfda/mpdp-routinesR.txt

    mpdp.obj <- mpdp(t(train_X),c(train_Y), nbmax=nbmax, nbbw=nbbw, Grid=c(1:nrow(train_X)),
      pcvpar=pcvpar,  kind.of.kernel="quadratic2")
    p.mpdp <- predict.mpdp(mpdp.obj, t(test_X))

    result[[1]][K,6] <- NA
    result[[2]][K,6] <- mean((c(test_Y) - p.mpdp)^2, na.rm=T)

    ##################################################
    ####################### NP #######################
    ### available at http://www.lsp.ups-tlse.fr/staph/npfda/npfda-routinesR.txt

    np.obj <- funopare.knn.gcv(c(train_Y), t(train_X), t(test_X), 0, nknot=nknot, c(0,1))
    p.np <- np.obj$Predicted.values

    result[[1]][K,7] <- NA
    result[[2]][K,7] <- mean((c(test_Y) - p.np)^2, na.rm=T)

    ##################################################
    #################### RDG #########################

    rdg.obj <- rdg(dat=dat[, 1:n1, drop=F])
    p.rdg <- pred.rdg(rdg.obj, newX=test_X)

    result[[1]][K,8] <- NA
    result[[2]][K,8] <- mean((p.rdg$yhat-c(test_Y))^2, na.rm=T)

  }

  return(result)

}

#######################################
#   1) Fertility rate
#   from Human Fertility Database
#   (www.humanfertility.org)
#######################################
load("fr.Rdata") # available from the github page, originally downloaded from www.humanfertility.org
ls()
data.fert <- t(dat[,,dimnames(dat)[[3]]==20])
result <- DATAanlys(dat0=data.fert, N=100, L=9, maxq=4, prtn=0.83, nknot=5, nbmax=5, nbbw=5, pcvpar=1/4, s=seq(0.00000001,0.00001,length.out=30), inisp=1)

### Table 3
MSPE <- rbind(
  round(apply(result$mspe, 2, mean)*(10^6), 2),
  round(apply(result$mspe, 2, median)*(10^6), 2),
  round(apply(result$mspe, 2, sd)*(10^6), 2)
)
colnames(MSPE) = c("MLR", "FLR", "FLiRTI", "SRPl", "SRPc", "MPDP", "NP", "RDG")
rownames(MSPE) <- c("Mean", "Median", "sd")
MSPE

### Figure 6
layout(matrix(c(1,1,2,2,3,3,4,4),nrow=1))
barplot(table(result$qhat[,1]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,80))
barplot(table(result$qhat[,4]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,100))
barplot(table(result$qhat[,5]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,80))
tr <- as.matrix(table(apply(result$mspe, 1, which.min)))
rownames(tr) <- c("MLR", "FLR", "FLiRTI", "SRP.L","SRP.C", "MPDP", "RDG")
barplot(t(tr), cex.axis=1, cex.names=1, ylim=c(0,30))



#######################################
#   2) Nitrogen oxides level
#   from the R package aire.zmvm
#######################################
library(aire.zmvm)
NOX <- get_station_data("HORARIOS", "NOX", 1986:2016)
NOX$DATE <- as.Date(NOX$date)
NOX <- separate(NOX, "DATE", c("Year", "Month", "Day"), sep = "-")
stlev <- unique(NOX$station_code)
NOX <- NOX[NOX$station_code==stlev[4],]

lev <- unique(NOX$date)
numobs <- rep(NA, length(lev))
for(i in 1:length(lev)){
  numobs[i] <- sum(NOX$date==as.numeric(lev[i]))
}
NOX <- NOX[! NOX$date %in% as.numeric(lev[which(numobs!=24)]),]
NOX <- NOX[NOX$Year==2016,]

dat <- matrix(NA, ncol=dim(NOX)[1]/24, nrow=24)
newlev <- unique(NOX$date)
for(i in 1:ncol(dat)){
  dat[,i] <- NOX$value[which(NOX$date==newlev[i])]
}
dat <- sqrt(dat)
# remove the curves which do not contain all the 24 hours observations
which(colSums(is.na(dat)) > 0)
dat <- dat[, -which(colSums(is.na(dat)) > 0)]

par(mfrow=c(1,1), mar=rep(3,4))
plot(dat[,1], col=2, type="l", ylim=range(dat, na.rm=T))
for(i in 2:dim(dat)[2]){
  lines(dat[,i], col=i)
}

result <- DATAanlys(dat0=dat, N=100, L=9, maxq=3,  prtn=0.65, nknot=5, nbmax=5, nbbw=5, pcvpar=1, s=seq(0.000001,0.001,length.out=30), inisp=1)

### Table 4
MSPE <- rbind(
  round(apply(result$mspe, 2, mean)*(10^2), 2),
  round(apply(result$mspe, 2, median)*(10^2), 2),
  round(apply(result$mspe, 2, sd)*(10^2), 2)
)
colnames(MSPE) = c("MLR", "FLR", "FLiRTI", "SRPl", "SRPc", "MPDP", "NP", "RDG")
rownames(MSPE) <- c("Mean", "Median", "sd")
MSPE

### Figure 9
layout(matrix(c(1,1,2,2,3,3,4,4),nrow=1))
barplot(table(result$qhat[,1]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,60))
barplot(table(result$qhat[,4]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,80))
barplot(table(result$qhat[,5]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,80))
tr <- as.matrix(table(apply(result$mspe, 1, which.min)))
rownames(tr) <- c("MLR", "FLR", "FLiRTI", "SRP.L","SRP.C", "MPDP", "RDG")
barplot(t(tr), cex.axis=1, cex.names=1, ylim=c(0,40))



#######################################
#   3) Disney stock volatility series
#   from Wharton Research Data Services
#   (wrds-web.wharton.upenn.edu/wrds/)
#######################################

result <- DATAanlys(dat0=dat, N=100, L=35, maxq=30,  prtn=0.5, nknot=20, nbmax=5, nbbw=5, pcvpar=1/4, s=seq(0.00005, 1, length.out=30), inisp=10)

### Table 5
MSPE <- rbind(
  round(apply(result$mspe, 2, mean), 2),
  round(apply(result$mspe, 2, median), 2),
  round(apply(result$mspe, 2, sd), 2)
)
colnames(MSPE) = c("MLR", "FLR", "FLiRTI", "SRPl", "SRPc", "MPDP", "NP", "RDG")
rownames(MSPE) <- c("Mean", "Median", "sd")
MSPE

### Figure 11
layout(matrix(c(1,1,2,2,3,3,4,4),nrow=1))
barplot(table(result$qhat[,1]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,60))
barplot(table(result$qhat[,4]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,50))
barplot(table(result$qhat[,5]), main="", cex.axis=1.5, cex.names=1.5, ylim=c(0,100))
tr <- as.matrix(table(apply(result$mspe, 1, which.min)))
rownames(tr) <- c("MLR","FLiRTI","SRP.L","SRP.C")
barplot(t(tr), cex.axis=1, cex.names=1, ylim=c(0,60))



#######################################
#   4) Montly sunspot data
#   from Base R datasets
#######################################
x <- sunspot.month
x <- sqrt(x)
x <- scale(x)

p <- 0.7
n <- length(x)
n1 <- floor(length(x)*p)
n2 <- n - n1
train.x <- x[1:n1]
test.x <- x[(n1+1):length(x)]

train.x <- as.ts(train.x)
test.x <- as.ts(test.x)
train.x.ts <- as.xts(train.x)
test.x.ts <- as.xts(test.x)

max.lag <- 150
lagged.train.x <- matrix(train.x.ts[-c(1:max.lag),1], ncol=1)
lagged.test.x <- matrix(test.x.ts[-c(1:max.lag),1], ncol=1)
for(i in 1:max.lag){
  lagged.train.x <- cbind(lagged.train.x, lag(train.x.ts, i)[-c(1:max.lag),1])
}
for(i in 1:max.lag){
  lagged.test.x <- cbind(lagged.test.x, lag(test.x.ts, i)[-c(1:max.lag),1])
}
colnames(lagged.train.x) <- c("lag0", paste0('lag',c(1:max.lag)))
colnames(lagged.test.x) <- c("lag0", paste0('lag',c(1:max.lag)))

x1 <- t(as.matrix(lagged.train.x[,dim(lagged.train.x)[2]:1]))
x2 <- t(as.matrix(lagged.test.x[,dim(lagged.test.x)[2]:1]))
dat <- cbind(x1, x2)

result <- list(qhat=matrix(NA, nrow=1, ncol=7), mspe=matrix(NA, nrow=1, ncol=7))
result <- lapply(result, function(x) {colnames(x) <- c("MLR", "FLR", "FLiRTI", "SRPl", "SRPc", "RDG", "OLS"); x})
p <- dim(dat)[1]
n <- dim(dat)[2]
n1 <- dim(x1)[2]
n2 <- dim(x1)[2]
train_X <- dat[c(1:(p-1)), 1:n1, drop=F]
test_X <- dat[c(1:(p-1)), (n1+1):n, drop=F]
train_Y <- dat[p, 1:n1, drop=F]
test_Y <- dat[p, (n1+1):n, drop=F]
##################################################
##################### MLR ########################
mlr.obj <- mlr(dat=dat[, 1:n1, drop=F], maxq=15)
p.mlr <- pred.mlr(mlr.obj, newX=test_X)
result[[1]][K,1] <- mlr.obj$qhat.mlr
result[[2]][K,1] <- mean((p.mlr$yhat-c(test_Y))^2, na.rm=T)
##################################################
##################### FLR ########################
flr.obj <- flr(dat=dat[, 1:n1, drop=F], L=35)
p.flr <- pred.flr(flr.obj, newX=test_X)
result[[1]][K,2] <- NA
result[[2]][K,2] <- mean((p.flr$yhat-c(test_Y))^2, na.rm=T)
##################################################
#################### FLiRTI ######################
### the result
### available at http://www-bcf.usc.edu/~gareth/research/flrti
d <- c(2,3,4)
s <- seq(0.000001, 0.1, length.out=30)
f.der2 <- flrti.cv(Y=c(train_Y), X=t(train_X), sigma=s, deriv=d[1], weight=0.1)
f.der3 <- flrti.cv(Y=c(train_Y), X=t(train_X), sigma=s, deriv=d[2], weight=0.1)
f.der4 <- flrti.cv(Y=c(train_Y), X=t(train_X), sigma=s, deriv=d[3], weight=0.1)
opt.d <- d[which.min(c(f.der2$error[which.min(f.der2$error)],
  f.der3$error[which.min(f.der3$error)],f.der4$error[which.min(f.der4$error)]))]
min.err <- c(which.min(f.der2$error), which.min(f.der3$error), which.min(f.der4$error))
opt.sigma <- s[min.err[opt.d-1]]
flrti.obj <- flrti(Y=c(train_Y), X=t(train_X), sigma=opt.sigma, deriv=opt.d, weight=0.1, plot=F)
p.flrti <- predict.flrti(flrti.obj, t(test_X))
result[[1]][K,3] <- NA
result[[2]][K,3] <- mean((p.flrti-c(test_Y))^2, na.rm=T)
##################################################
###################### SRPl ######################
### use package "srp"
srpl.obj <- srp.l(dat[, 1:n1, drop=F], maxq=15, plot=F)
p.srpl <-  predict(srpl.obj, x=test_X)
result[[1]][K,4] <- srpl.obj$qhat
result[[2]][K,4] <- mean((p.srpl$yhat - c(test_Y))^2, na.rm=T)
##################################################
###################### SRPc ######################
### use package "srp"
srpc.obj <- srp.c(dat[, 1:n1, drop=F], L=35, maxq=15, inisp=1, plot=F)
p.srpc <-  predict(srpc.obj, x=test_X)
result[[1]][K,5] <- srpc.obj$qhat
result[[2]][K,5] <- mean((p.srpc$yhat - c(test_Y))^2, na.rm=T)
##################################################
#################### RDG #########################
rdg.obj <- rdg(dat=dat[, 1:n1, drop=F])
p.rdg <- pred.rdg(rdg.obj, newX=test_X)
result[[1]][K,6] <- NA
result[[2]][K,6] <- mean((p.rdg$yhat-c(test_Y))^2, na.rm=T)
##################################################
#################### OLS #########################
ols.obj <- lm(lag0 ~ .-1, data=lagged.train.x)
lagged.test.x0 <- lagged.test.x[, -1, drop=F]
p.ols <- predict(ols.obj, newdata = lagged.test.x0)
result[[1]][K,7] <- NA
result[[2]][K,7] <- mean((p.ols-lagged.test.x$lag0)^2, na.rm=T)

### Figure 13
par(mfrow=c(1,3), mar=c(4.5,4,2,2))
plot(dat,  type="l", xlab="Year", ylab="")
pacf(as.vector(dat), lag.max=150, xlab="Lag", main="")
acf(as.vector(dat), lag.max=150, xlab="Lag", main="")

### Table 6
round(result$mspe, 4)*100
