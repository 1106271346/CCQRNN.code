library(qrnn)
library("nleqslv") 
library("quantreg")
###-------- DAU --------###
Drawy <- function(X, y, delta, U, L, B){
    # Only allows right censoring...
    cen <- which(delta == 0)
    if(length(cen)){
    for (j in 1:length(cen)) {
	yj <- B %*% X[cen[j]] 
	if (min(yj) < L[cen[j]]) {
           lo <- which(yj < L[cen[j]])
           if (length(lo) != 1) 				
              y[cen[j]] <- yj[sample(lo, 1)]
	   else
              y[cen[j]] <- yj[lo]
  	   }
	}
    }
    y
}
Draweta <- function(X, y, delta, Z, B,grid = 1:99/100, gam){
    n <- length(y)
    wts <- plogis(c(Z %*% gam))
    fit <- X %*% t(B)
    ## Computes the quantile level  where X falls based on initial estimator of beta
    indvec <- sapply((1:n), function(u) round(rank(c(y[u],fit[u,]))[1]))
    egrid <- c(grid,1) 
    qsurv <- (1 - egrid[indvec])
    tmp <- wts*qsurv
    hatp <- delta+(1-delta)*tmp/(1-wts+tmp)
    hatp <- ifelse(hatp>=0, hatp, 1e-6)
    rbinom(n, 1, hatp)
}

Drawgam <- function(Z, eta, bootstrap){
    n = length(eta)
    if (bootstrap) { 
	s <- sample(1:n, n, replace=TRUE)
	gam <- as.vector(glm(eta[s]~ 0 + Z[s,], binomial)$coef)
    }
    else 
	gam = as.vector(glm(eta ~ Z - 1,binomial)$coef)
    gam
}

DrawB <- function(yh, X,grid,bootstrap){
    n = length(yh)
    if (bootstrap) { 
	s <- sample(1:n, n, replace=TRUE)
	B <- t(rq(yh[s] ~ 0+X[s], tau=grid)$coef)
    }  
    else
	B <- t(rq(yh ~ 0+X, tau = grid)$coef)
    B
}

cqr.fit.DA <- function(X, Z, y, delta, link, grid = 1:99/100, B = NULL, gam=NULL, 
	   taus, L=NULL, U=NULL, h = NULL, bootstrap = TRUE, maxit = 50, Large = 1e4){
 yy=Y[train];XX=X[train];ZZ=Z[train,];delta1=delta[train]   
 n <- length(yy)
    if(is.null(B)){
	uncen <- which(delta1 == 0)
	f <- rq(yy~ XX - 1, tau=grid, subset = uncen)$coef  
      B=t(f) 
    }
    if(is.null(gam))
	gam = glm(delta1 ~ ZZ - 1, binomial(link = link))$coef
    if(is.null(L)) L = rep(-Large,n) 
    if(is.null(U)) {
      cen = which(delta1==0) 
      U = 0*L + Large 
      U[cen] = yy[cen]
    }
    G <- matrix(0, length(taus), maxit)
    H <- matrix(0, ncol(ZZ), maxit)
    h <- 0
    MSE1=MSE2=MSE3=MSE4=MSE5=MAE1=MAE2=MAE3=MAE4=MAE5=MSE_com=MAE_com=c()
    while (h < maxit){
	h <- h + 1
	eta <- Draweta(XX, yy, delta1, ZZ, B, grid, gam)
	gam <-Drawgam (ZZ, eta,bootstrap)
	H[,h] <- gam
	uncured <- which(eta==1)
	yh <- Drawy(XX[uncured], yy[uncured], delta1[uncured], U[uncured], L[uncured], B)
	B <- DrawB(yh, XX[uncured],grid,bootstrap)
	if(sum(B^2) < 10000 * length(grid)) 
	    G[,h] <- rq(yh ~ XX[uncured] - 1, tau = taus)$coef
	else G[,h] <- matrix(NA, ncol(XX),length(taus))
    XXX=X[-train];   
    yhat1=G[1,h]*(XXX[which(delta[-train]==1)]);yhat2=G[2,h]*(XXX[which(delta[-train]==1)]);
    yhat3=G[3,h]*(XXX[which(delta[-train]==1)]);
    yhat4=G[4,h]*(XXX[which(delta[-train]==1)]);yhat5=G[5,h]*(XXX[which(delta[-train]==1)])

    MSE1[h]=mean((ytrue1-yhat1)^2);MAE1[h]=mean(abs(ytrue1-yhat1))
    MSE2[h]=mean((ytrue2-yhat2)^2);MAE2[h]=mean(abs(ytrue2-yhat2))
    MSE3[h]=mean((ytrue3-yhat3)^2);MAE3[h]=mean(abs(ytrue3-yhat3))
    MSE4[h]=mean((ytrue4-yhat4)^2);MAE4[h]=mean(abs(ytrue4-yhat4))
    MSE5[h]=mean((ytrue5-yhat5)^2);MAE5[h]=mean(abs(ytrue5-yhat5))
    
    ytrue_com=c(ytrue1,ytrue2,ytrue3,ytrue4,ytrue5)
    yhat_com=c(yhat1,yhat2,yhat3,yhat4,yhat5)
    MSE_com[h]=mean((ytrue_com-yhat_com)^2);MAE_com[h]=mean(abs(ytrue_com-yhat_com))
  }
  list(beta = apply(G, 1, mean), gamma = apply(H, 1, mean),MSE=c(mean(MSE_com),mean(MSE1),mean(MSE2),mean(MSE3),mean(MSE4),mean(MSE5))
      ,MAE=c(mean(MAE_com),mean(MAE1),mean(MAE2),mean(MAE3),mean(MAE4),mean(MAE5)))
} 

###-------- IMP --------###
impute <-function(X, Z, y, D,g, taus, link, nimp = 5, grid = NULL, h = NULL){
    n <- length(y[train])
    lam <- rep(1,n)
    Tau <- max(y[train] * D[train])
    B <- NWweight(Z[train,], h, lam)
    G <- gfit(Tau, y[train], D[train], Z[train,], B, g$coef, lam)
    pZg <-plogis(Z[train,] %*% G$g)
    S= diag(G$S)
    weights = pZg

    if(is.null(grid))
	grid <- seq(0,1,by=min(0.01,1/(2*length(y)^.7)))
    if (any(weights < 0)) stop("negative weights")
    KMW <- function(a,b) 
	ifelse( (b-1+a > 0 & a >0), (b-1+a)/a, 1) 
    p <- 1
    n <- length(y[train])
    m <- length(grid)
    Tau <- max(y[train] * D[train])
    tmp <- weights * S
    hatp <- D[train] + (1 - D[train]) * tmp/ifelse(1 - weights + tmp == 0, 1, 1 - weights + tmp)
    hatp <- pmax(hatp, 1e-6)
    hbet <- matrix(0,p,m)
    xbar <- mean(X[train])
    for(i in 1:nimp){ 
	eta_train = rbinom(n[train],1,hatp)
	omeF <- outer(S,grid, KMW)
	omeF[which(D[train]==1),] <- 1
	uncr <- which(eta_train == 1)
	kmweights <- omeF[uncr,]; 
	y.uncr <- y[uncr]; D.uncr <- D[uncr]; X.uncr <- X[uncr];lam.uncr <- lam[uncr]
	# censored objects among uncured
	ncen <- sum(D.uncr == 0)
	y.aug <- rep(1e6, ncen) ## Points augmented at infinity as in Portnoy's CRQ
	X.aug <-  X.uncr[D.uncr==0] 
	kmweights <- lam.uncr*kmweights
	kmweights.aug <- lam.uncr[D.uncr==0]*(1-kmweights[D.uncr ==0,])
	ynew <- c(y.uncr,y.aug)  
	Xnew <- c(X.uncr,X.aug)
	wts <- rbind(kmweights,kmweights.aug) 
	hbet <- hbet + sapply((1:m), function(j) rq(ynew ~ 0+Xnew,tau=grid[j],weights = wts[,j])$coef)
   }
   beta <- hbet/nimp
   Qhat <- t(xbar) %*% beta
   beta <- rbind(grid,beta, Qhat)
   class(beta) <- "cqr"
    est.imp.beta= coef(beta, taus)
    XXX=X[-train];
    yhat1= est.imp.beta[1]*XXX[which(D[-train]==1)];
    yhat2= est.imp.beta[2]*XXX[which(D[-train]==1)];
    yhat3= est.imp.beta[3]*XXX[which(D[-train]==1)];
    yhat4= est.imp.beta[4]*XXX[which(D[-train]==1)];
    yhat5= est.imp.beta[5]*XXX[which(D[-train]==1)];
    ytrue_com=c(ytrue1,ytrue2, ytrue3,ytrue4,ytrue5)
    yhat_com=c(yhat1,yhat2,yhat3,yhat4,yhat5)
    MSE1=mean((ytrue1-yhat1)^2);MSE2=mean((ytrue2-yhat2)^2);MSE3=mean((ytrue3-yhat3)^2);
    MSE4=mean((ytrue4-yhat4)^2);MSE5=mean((ytrue5-yhat5)^2);
    MAE1=mean(abs(ytrue1-yhat1));MAE2=mean(abs(ytrue2-yhat2));MAE3=mean(abs(ytrue3-yhat3));
    MAE4=mean(abs(ytrue4-yhat4));MAE5=mean(abs(ytrue5-yhat5))
    MSE_com=mean((ytrue_com-yhat_com)^2); MAE_com=mean(abs(ytrue_com-yhat_com));
return(list(beta=est.imp.beta,MSE=c(MSE_com,MSE1,MSE2,MSE3,MSE4,MSE5),MAE=c(MAE_com,MAE1,MAE2,MAE3,MAE4,MAE5),gamma=G$g))
}
coef.cqr <- function (object, taus = 1:4/5, ...) {
    if (min(taus) < 0 || max(taus) > 1) 
        stop("taus out of range [0,1]")
    taus <- sort(taus)
    S <- object
    r <- S[1, ]
    r <- c(r[1], r)
    if (is.unsorted(r)) 
        r <- r[-length(r)]
    B <- S[-1, , drop = FALSE]
    B <- t(cbind(B, B[, ncol(B), drop = FALSE]))
    ts <- taus[taus > min(r) & taus < max(r)]
    bin <- findInterval(ts, r)
    wgt <- (ts - r[bin])/(r[bin + 1] - r[bin])
    binlag <- bin - 1
    binlag[binlag == 0] <- 1
    coef <- t(wgt * B[bin, , drop = FALSE] + 
	 (1 - wgt) * B[binlag, , drop = FALSE])
    nna <- length(taus) - length(ts)
    if (nna > 0) 
        coef <- cbind(coef, matrix(NA, nrow(coef), nna))
    taulabs <- paste("tau=", format(round(taus, 3)))
    dimnames(coef)[[2]] <- taulabs
    coef[-nrow(coef), ]
}

biweight <- function(x) (15/16) * (abs(x) <= 1) * (1 - x^2)^2
NWweight <- function(Z, h, Lam, K = biweight){
	n <- NROW(Z)
	p <- NCOL(Z)
	B <- apply(Z, 2, function(x) outer(x,x,"-"))
	B <- array(B, c(n,n,p))
	for(i in 1:p) B[,,i] <- K(B[,,i]/h[i])
	B <- apply(B, 1:2, prod)
	B/apply(B,1,sum)
}
LocalKM <- function(y, D, B, tau, omega = NULL, eta = NULL){
	n <- length(y)
	o <- order(y)
	y <- y[o]
	D <- D[o]
	B <- B[o,o]
	if(length(eta)) B <- eta[o] * B
	N <- B %*% (D * outer(y, y, '<='))
	dN <- N - cbind(0, N[,-n])
	if(length(omega)) B <- omega[o,o] * B
	R <- B %*% outer(y, y, '>=')
	S <- 1 - dN/R
	S[is.nan(S)] <- 1
	S <- apply(S, 1, cumprod)
	S[y >= tau,] <- 0
	return(S[order(o),order(o)])
}
Omega <- function(D, Z, g, S){
	p <- plogis(c(Z %*% g))
	pS <- p * t(S)
	t(D + (1 - D) * pS/(1 - p + pS))
}
gfit <- function(tau, y, D, Z, B, g0, Lam, cntl = list(btol = 0.01)){
    n <- nrow(Z)
    p <- ncol(Z)
    omega <- matrix(1,n,n)
    S <- LocalKM(y, D, B, tau, omega)
    EEg <- function(g){
	p <- plogis(c(Z %*% g))
	R <- p * (1 - diag(S))
	R <- Z * Lam * (1 - p) * (D - R)/(1 - R)
	apply(R,2,sum)
    }
    gap <- 1
    while(gap > 1e-2){
	omega <- Omega(D, Z, g0, S)
	S <- LocalKM(y, D, B, tau, omega)
	fitg <- nleqslv(g0, EEg, control = cntl)
	eflag <- fitg$termcd
	if(eflag != 1) break
	gap <- max(abs(g0 - fitg$x))
	g0 <- fitg$x
    }
    omega <- Omega(D, Z, g0, S)
    S <- LocalKM(y, D, B, tau, omega)
    list(g = g0, S = S, eflag = eflag)
}
###-------- CCQRNN-------- ###
Drawgamc <- function(Z, eta, bootstrap){
  n = length(eta)
  if (bootstrap) { 
    s <- sample(1:n, n, replace=TRUE)
    gam <- as.vector(glm(eta[s]~ 0 + Z[s,], binomial(link = logit))$coef)
  }
  else 
    gam = as.vector(glm(eta ~ Z - 1,binomial)$coef)
  gam
}

DrawB.parmsc <- function(yh, X,grid,bootstrap){
  n = length(yh)
  if (bootstrap) { 
    s <- sample(1:n, n, replace=TRUE)
    
    
    B.parms= mcqrnn.fit(as.matrix(X[s]), as.matrix(yh[s]), n.hidden=5, w=NULL, tau=grid,
                        iter.max=5000, n.trials=1,
                        lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5),
                        eps.seq=2^seq(-8, -32, by=-4), Th=sigmoid,
                        Th.prime=sigmoid.prime, penalty=0, n.errors.max=10,
                        trace=TRUE)
      }  
  B.parms
}

Drawetac<- function(X, y, delta, Z, B.parms,grid = 1:99/100, gam){
  n <- length(y)
  wts <- plogis(c(Z %*% gam))
  
  fit=mcqrnn.predict(as.matrix(X), B.parms)  
  indvec <- sapply((1:n), function(u) round(rank(c(y[u],fit[u,]))[1]))
  egrid <- c(grid,1) 
  qsurv <- (1 - egrid[indvec])
  tmp <- wts*qsurv
  hatp <- delta+(1-delta)*tmp/(1-wts+tmp)
  hatp <- ifelse(hatp>=0, hatp, 1e-6)
  rbinom(n, 1, hatp)
}
Drawyc <- function(X, y, delta, U, L, B.parms){
  cen <- which(delta == 0)
  if(length(cen)){
    for (j in 1:length(cen)) {
      yj <- mcqrnn.predict(as.matrix(X[cen[j]]), B.parms)
      if (min(yj) < L[cen[j]]) {
        lo <- which(yj < L[cen[j]])
        if (length(lo) != 1) 				
          y[cen[j]] <- yj[sample(lo, 1)]
        else
          y[cen[j]] <- yj[lo]
      }
    }
  }
  y
}

cqr.fit.DAc <- function(X, Z, y, delta, link, grid = 1:99/100, B = NULL, gam=NULL, 
                       taus, L=NULL, U=NULL, h = NULL, bootstrap = TRUE, maxit =50, Large = 1e4,I11,lambda11){
  yy=Y[train];XX=X[train];ZZ=Z[train,];delta1=delta[train]
  n <- length(yy)
  if(is.null(B)){
    uncen <- which(delta1 == 0)
    
    B.parms= mcqrnn.fit(as.matrix(XX[uncen]), as.matrix(yy[uncen]), n.hidden=5, w=NULL, tau=grid,
                        iter.max=5000, n.trials=1,
                        lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5),
                        eps.seq=2^seq(-8, -32, by=-4), Th=sigmoid,
                        Th.prime=sigmoid.prime, penalty=0, n.errors.max=10,
                        trace=TRUE)
  }  
  
  if(is.null(gam))
    gam = glm(delta1 ~ ZZ-1 , binomial(link = logit))$coef
  if(is.null(L)) L = rep(-Large,n) 
  if(is.null(U)) {
    cen = which(delta1==0) 
    U = 0*L + Large 
    U[cen] = yy[cen]
  }
  G <- array(0, c(ncol(X), length(taus), maxit)) 
  H <- matrix(0, ncol(ZZ), maxit)
  h <- 0
  ETA=matrix(0,maxit,n)
  YY=matrix(0,maxit,length(taus))
  B.final=vector("list",maxit)
  MSE1=MSE2=MSE3=MSE4=MSE5=MAE1=MAE2=MAE3=MAE4=MAE5=MSE_com=MAE_com=c()
  while (h < maxit){
    h <- h + 1
    eta11 <- Drawetac(XX, yy, delta1, ZZ, B.parms, grid, gam)
    gam <-Drawgamc (ZZ, eta11,bootstrap)
    H[,h] <- gam
    uncured <- which(eta11==1)
    yh <- Drawyc(XX[uncured], yy[uncured], delta1[uncured], U[uncured], L[uncured], B.parms)
    B.parms <- DrawB.parmsc(yh, XX[uncured],grid,bootstrap)
    ETA[h,]=eta11
    B.final[[h]]=mcqrnn.fit(as.matrix(XX[uncured]), as.matrix(yh), n.hidden=I11, w=NULL, tau=taus,
                            iter.max=5000, n.trials=1,
                            lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5),
                            eps.seq=2^seq(-8, -32, by=-4), Th=sigmoid,
                            Th.prime=sigmoid.prime, penalty=lambda11, n.errors.max=10,
                            trace=TRUE)
    
    XXX=X[-train]
    G=mcqrnn.predict(as.matrix(XXX[delta[-train]==1]), B.final[[h]])
    
    YY[h,]=apply(G,2,mean)    
    MSE1[h]=mean((ytrue1-G[,1])^2); MAE1[h]=mean(abs(ytrue1-G[,1]))
    MSE2[h]=mean((ytrue2-G[,2])^2); MAE2[h]=mean(abs(ytrue2-G[,2]))
    MSE3[h]=mean((ytrue3-G[,3])^2); MAE3[h]=mean(abs(ytrue3-G[,3]))
    MSE4[h]=mean((ytrue4-G[,4])^2); MAE4[h]=mean(abs(ytrue4-G[,4]))
    MSE5[h]=mean((ytrue5-G[,5])^2); MAE5[h]=mean(abs(ytrue5-G[,5]))
    
    ytrue_com=c(ytrue1,ytrue2,ytrue3,ytrue4,ytrue5)
    MSE_com[h]=mean((ytrue_com-c(G))^2); MAE_com[h]=mean(abs(ytrue_com-c(G)))
    print(MSE_com)
}
  list(gamma = apply(H, 1, mean),MSE=c(mean(MSE_com),mean(MSE1),mean(MSE2),mean(MSE3),mean(MSE4),mean(MSE5))
       ,MAE=c(mean(MAE_com),mean(MAE1),mean(MAE2),mean(MAE3),mean(MAE4),mean(MAE5)))
} 

cqr.fit.DA.pre <- function(X, y, delta, taus, I11,lambda11,B=NULL){
  n <- length(y)
  if(is.null(B)){
    uncen <- which(delta == 0)
    
    B.parms= mcqrnn.fit(as.matrix(X[uncen]), as.matrix(y[uncen]), n.hidden=I11, w=NULL, tau=taus,
                        iter.max=5000, n.trials=1,
                        lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5),
                        eps.seq=2^seq(-8, -32, by=-4), Th=sigmoid,
                        Th.prime=sigmoid.prime, penalty=lambda11, n.errors.max=10,
                        trace=TRUE)
  }
  
  list(parms= B.parms)
}
par.choice=function(tau,y,Z,n){
Kv=5;n1=seq(0,n,n/Kv)
parms=p=vector("list",length(tau))
R=8;M=8;K=length(tau)
I=seq(1,8,length=R)
lambda=seq(0.01,0.08,length=M)
MSE_CV=matrix(0,R,M); MSES=c()
pp.final=matrix(0,Kv,K*n)
y1=bet0*Z+qerror1
y2=bet0*Z+qerror2
y3=bet0*Z+qerror3
y4=bet0*Z+qerror4
y5=bet0*Z+qerror5
yy=cbind(y1,y2,y3,y4,y5)
  for(r in 1:R){
    for(m in 1:M){
      I1=I[r]
      lambda1=lambda[m]
      for(i in 1:Kv){
        y_train=Y[-c((n1[i]+1):n1[i+1])]
        X_train=Z[-c((n1[i]+1):n1[i+1])]
        delta_train=delta[-c((n1[i]+1):n1[i+1])]
        parms=cqr.fit.DA.pre(X=X_train, y=y_train, delta=delta_train, taus=taus,I11=I1,lambda11=lambda1,B=NULL)$parms
        pp=mcqrnn.predict(as.matrix(Z[c((n1[i]+1):n1[i+1])]), parms)
        pp.final=c(pp)
        yyy=yy[c((n1[i]+1):n1[i+1]),]
        #predict=apply(pp.final,2,mean)
        MSES[i]=mean((pp.final-c(yyy))^2)
      }
      MSE_CV[r,m]=mean(MSES)
    }
  } 
min_value=min(MSE_CV)
lo=which(MSE_CV == min_value, arr.ind = TRUE)
I11=I[lo[1]];lambda11=lambda[lo[2]]
return(list(I11=I11,lambda11=lambda11))
}

###-------- CQRNN-------- ###
DArq=function(X, y, delta, link, grid = 1:99/100, B = NULL, 
                       taus, L=NULL, U=NULL, h = NULL, bootstrap = TRUE, maxit =50, Large = 1e4,I11,lambda11){
  yy=y[train];XX=X[train];delta1=delta[train]
  n <- length(yy)
  if(is.null(B)){
    uncen <- which(delta1 == 0)
    
    B.parms= mcqrnn.fit(as.matrix(XX[uncen]), as.matrix(yy[uncen]), n.hidden=5, w=NULL, tau=grid,
                        iter.max=5000, n.trials=1,
                        lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5),
                        eps.seq=2^seq(-8, -32, by=-4), Th=sigmoid,
                        Th.prime=sigmoid.prime, penalty=0, n.errors.max=10,
                        trace=TRUE)
  }  
  
  if(is.null(L)) L = rep(-Large,n) 
  if(is.null(U)) {
    cen = which(delta1==0) 
    U = 0*L + Large 
    U[cen] = yy[cen]
  }
  G <- array(0, c(ncol(X), length(taus), maxit)) 
  h <- 0
  YY=matrix(0,maxit,length(taus))
  B.final=vector("list",maxit)
  MSE1=MSE2=MSE3=MSE4=MSE5=MAE1=MAE2=MAE3=MAE4=MAE5=MSE_com=MAE_com=c()
  while (h < maxit){
    h <- h + 1
    yh <- Drawyc(XX, yy, delta1, U, L, B.parms)
    B.parms <- DrawB.parmsc(yh, XX,grid,bootstrap)
    B.final=mcqrnn.fit(as.matrix(XX), as.matrix(yh), n.hidden=I11, w=NULL, tau=taus,
                            iter.max=5000, n.trials=1,
                            lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5),
                            eps.seq=2^seq(-8, -32, by=-4), Th=sigmoid,
                            Th.prime=sigmoid.prime, penalty=lambda11, n.errors.max=10,
                            trace=TRUE)
    
    
    XXX=X[-train]
    G=mcqrnn.predict(as.matrix(XXX[delta[-train]==1]), B.final)
    
    YY[h,]=apply(G,2,mean)
    
    MSE1[h]=mean((ytrue1-G[,1])^2); MAE1[h]=mean(abs(ytrue1-G[,1]))
    MSE2[h]=mean((ytrue2-G[,2])^2); MAE2[h]=mean(abs(ytrue2-G[,2]))
    MSE3[h]=mean((ytrue3-G[,3])^2); MAE3[h]=mean(abs(ytrue3-G[,3]))
    MSE4[h]=mean((ytrue4-G[,4])^2); MAE4[h]=mean(abs(ytrue4-G[,4]))
    MSE5[h]=mean((ytrue5-G[,5])^2); MAE5[h]=mean(abs(ytrue5-G[,5]))
    
    ytrue_com=c(ytrue1,ytrue2,ytrue3,ytrue4,ytrue5)
    MSE_com[h]=mean((ytrue_com-c(G))^2); MAE_com[h]=mean(abs(ytrue_com-c(G)))
    print(MSE_com)
}
  list(MSE=c(mean(MSE_com),mean(MSE1),mean(MSE2),mean(MSE3),mean(MSE4),mean(MSE5))
       ,MAE=c(mean(MAE_com),mean(MAE1),mean(MAE2),mean(MAE3),mean(MAE4),mean(MAE5)))
} 

###-------- main code-------- ### 
n=400      # sample size
Large=1e5;
taugrid = seq(0.01, 0.99, by=0.01);
logitf<-function(u){exp(u)/(1+exp(u))}
# true value of the parameters
gam0=c(1,-1);p2=length(gam0)
bet0=c(1);p1=length(bet0);
colZ=c(4,5)
taus=c(0.3,0.4,0.5,0.6,0.7);K=length(taus)
sigma=0.2;
qerror1=sigma*qnorm(taus[1],0,0.5)
qerror2=sigma*qnorm(taus[2],0,0.5)
qerror3=sigma*qnorm(taus[3],0,0.5)
qerror4=sigma*qnorm(taus[4],0,0.5)
qerror5=sigma*qnorm(taus[5],0,0.5)

W=cbind(1,rnorm(n,0,0.25));
Z=W[,2];

gamW=c(W%*%gam0);
peta=logitf(gamW)
eta=rbinom(n,1,peta)

error=sigma*rnorm(n,0,0.5)
TT=bet0*Z+error;
T=(1-eta)*Large+eta*TT;

L=3;U=runif(n,0,L);
C=U*(Z<0.5)+(U+1)*(Z>=0.5)

Y=pmin(T,C)
delta=as.numeric(T<=C)
dat=cbind(Y,delta,eta,Z)
id=c(1:n)
train=sample(id,n*0.8,replace=FALSE)

XXX=Z[-train]
ytrue1=bet0*(XXX[which(delta[-train]==1)])+qerror1
ytrue2=bet0*(XXX[which(delta[-train]==1)])+qerror2
ytrue3=bet0*(XXX[which(delta[-train]==1)])+qerror3
ytrue4=bet0*(XXX[which(delta[-train]==1)])+qerror4
ytrue5=bet0*(XXX[which(delta[-train]==1)])+qerror5


DA=cqr.fit.DA(X=Z, Z=W, y=Y, delta, link="logit", grid = 1:99/100, B = NULL, gam=NULL, 
            taus, L=NULL, U=NULL, h = NULL, bootstrap = TRUE, maxit = 50, Large = 1e5)

g2=(glm(delta~Z,binomial))
IMP =impute(X=Z, Z=W,y=Y, D=delta, g=g2, taus, link="logit", nimp = 5,grid = 1:99/100, h = rep(0.6027438,p2))

par=par.choice(taus,y=Y,Z=Z,n);I11=par$I11;lambda11=par$lambda11
NN=cqr.fit.DAc(X=Z, Z=W, y=Y, delta, link, grid = 1:99/100, B = NULL, gam=NULL, 
                    taus, L=NULL, U=NULL, h = NULL, bootstrap = TRUE, maxit =50, Large = 1e4,I11=I11,lambda11=lambda11)



nocure=DArq(X=Z, y=Y, delta, link, grid = 1:99/100, B = NULL, 
                    taus, L=NULL, U=NULL, h = NULL, bootstrap = TRUE, maxit =50, Large = 1e4,I11=I11,lambda11=lambda11)


hat.gamma=c(IMP$gamma,DA$gamma,NN$gamma)
hat.beta=c(IMP$beta,DA$beta)
NN_eta=1-mean(logitf(W[-train,]%*%NN$gamma));
DA_eta=1-mean(logitf(W[-train,]%*%DA$gamma));
IMP_eta=1-mean(logitf(W[-train,]%*%IMP$gamma))
MSE_hat=c(IMP$MSE,DA$MSE,NN$MSE,nocure$MSE)
MAE_hat=c(IMP$MAE,DA$MAE,NN$MAE,nocure$MAE)


