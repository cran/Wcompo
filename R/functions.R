
# data processing
# id: A numeric vector denoting the subject
# time: Time to death or censoring
# status:0=censoring; 1=death; 2, 3,...=different types of recurrent event
# Z: Design matrix
# w: Weight vector for event types 1, 2, ..., K
# ep: an infinitesimal tolerance threshold for error in event times
data.process.w=function(id,time,status,Z,w,ep)
{

  uid <- unique(id)
  n   <- length(uid)
  t   <- unique(time)
  t   <- t[order(t)]
  m   <- length(t)
  Z0  <- as.matrix(Z[status==0|status==1,])

  E         <- list()
  length(E) <- n
  X     <- rep(NA,n)
  delta <- rep(NA,n)
  dN    <- matrix(0,n,m)
  dNC   <- matrix(0,n,m)
  Y     <- matrix(0,n,m)
  YC    <- matrix(NA,n,m)
  for (i in 1:n)
  {


      E[[i]] <- cbind(time[id==uid[i]],
                      status[id==uid[i]])

      last     <- which(E[[i]][,2] == 0 | E[[i]][,2] == 1)
      X[i]     <- E[[i]][last,1]
      delta[i] <- as.numeric(E[[i]][last,2])

      #make E a matrix
      E[[i]] <- matrix(E[[i]], length(E[[i]])/2, 2)

      Y[i,X[i] >= t-ep] <- 1

      if (delta[i]==1){YC[i,]=1}else{YC[i,]=Y[i,]}

      K=length(w)


      for (j in 1:m)
      {
          for (k in 1:K)
          {
              dN[i,j] = dN[i,j] + w[k] * sum(abs(E[[i]][,1][E[[i]][,2]==k]-t[j]) < ep)
          }
      }
  }

  return(list(t=t,uid=uid,Z=Z0,E=E,X=X,delta=delta,dN=dN,dNC=dNC,Y=Y,YC=YC))
}




#' Fit a proportional means regression model for weighted composite endpoint
#' of recurrent event and death
#'
#' @description Fit a semiparametric proportional means regression model for the weighted
#' composite endpoint of recurrent event and death (Mao and Lin, 2016).
#' (Jared D. Huling (ORCID: 0000-0003-0670-4845) contributed to the optimization of this code.)
#'
#' @param id A vector of unique patient identifiers.
#' @param time A vector of event times.
#' @param status A vector of event type labels. 0: censoring; 1: death;
#'  2, 3,..., \eqn{K}: different types of (possibly recurrent) nonfatal event.
#' @param Z Covariate matrix (must be time-constant).
#' @param w A \eqn{K}-vector of weights assigned to event types 1 (death), 2, ..., \eqn{K}
#' (nonfatal events); If \code{NULL}, an unweighted endpoint is modeled
#' (i.e., with \code{w=c(1, 1, ..., 1)}).
#' @param ep Convergence threshold for the Newton-Raphson algorithm.
#' @return An object of class \code{CompoML} with the following components.
#' \code{beta}: a vector of estimated regression coefficients (log-mean ratios);
#' \code{var}: estimated covariance matrix for \code{beta};
#' \code{t}: unique event times;
#' \code{y}: estimated baseline mean function (of \code{t}).
#' @seealso \code{\link{plot.CompoML}}, \code{\link{print.CompoML}}
#' @export
#' @keywords CompoML
#' @references Mao, L. and Lin, D. Y. (2016). Semiparametric regression for the weighted
#' composite endpoint of recurrent and terminal events. Biostatistics, 17, 390-403.
#' @examples
#' \donttest{
#' ## load package and data
#' library(Wcompo)
#' head(hfmock)
#' ## fit a weighted PM (w_D=2, w_1=1)
#' obj <- CompoML(hfmock$id,hfmock$time,hfmock$status,hfmock[,c("Training","HF.etiology")],
#'                w=c(2,1))
#' ## print out the result
#' obj
#'
#' oldpar <- par(mfrow = par("mfrow"))
#' par(mfrow=c(1,2))
#' ## plot the estimated mean function for
#' ## non-ischemic patients by treatment
#' plot(obj,c(1,0),ylim=c(0,1.5),xlim=c(0,50),
#'      main="Non-ischemic",
#'      xlab="Time (months)",cex.main=1.2,lwd=2)
#' plot(obj,c(0,0),add=TRUE,cex.main=1.2,lwd=2,lty=2)
#' legend("topleft",lty=1:2,lwd=2,c("Exercise training","Usual care"))
#'
#'
#' ## plot the estimated mean function for
#' ## ischemic patients by treatment
#' plot(obj,c(1,1),ylim=c(0,1.5),xlim=c(0,50),
#'      main="Ischemic",
#'      xlab="Time (months)",cex.main=1.2,lwd=2)
#' plot(obj,c(0,1),add=TRUE,cex.main=1.2,lwd=2,lty=2)
#' legend("topleft",lty=1:2,lwd=2,c("Exercise training","Usual care"))
#' par(oldpar)
#' }
CompoML <- function(id,time,status,Z,w=NULL,ep=1e-4)
{

  # if w is NULL, then use the unweighted composite
  if (is.null(w)){
    K=length(unique(status))-1
    w=rep(1,K)
  }

  # pre-processing data
  Z=as.matrix(Z)
    pced.data <- data.process.w(id,time,status,Z,ep=ep,w=w)

    t     <- pced.data$t
    uid   <- pced.data$uid
    Z     <- as.matrix(pced.data$Z)
    X     <- pced.data$X
    delta <- pced.data$delta
    dN    <- pced.data$dN
    dNL   <- pced.data$dNC
    Y     <- pced.data$Y
    YL    <- pced.data$YC
    n     <- nrow(dN)
    m     <- ncol(dN)



    # fit Cox model for the censoring distribution
    # to get the IPCW weights
    coxph <- coxph.LTFU(t=t,Z=Z,X=X,delta=delta,dNL=dNL,Y=Y,YL=YL,ep=ep)
    dML   <- coxph$dML
    dL    <- coxph$dL
    W     <- coxph$W
    gamma <- coxph$gamma

    p   <- ncol(Z)

    # fit the proportional means model for the weighted composite
    # using post-processed data and IPCW weights
    fit <- NRbeta(W    = W,
                  Z    = Z,
                  dN   = dN,
                  init = rep(0, p),
                  err  = ep)

    beta <- fit$beta
    U    <- fit$U
    A    <- fit$A
    iter <- fit$i

    # variance estimation

    mart <- event.mart(beta=beta,W=W,Z=Z,dN=dN)
    du   <- mart$du
    dM   <- mart$dM
    S0   <- mart$S0
    S1   <- mart$S1
    Zb   <- mart$Zb

    Om    <- Omega(gamma=gamma,Y=Y,Z=Z,dNL=dNL)
    Omega <- Om$Omega
    R0    <- Om$R0
    R1    <- Om$R1

    bgq   <- BGQ(gamma=gamma,Z=Z,R0=R0,R1=R1,dL=dL,dM=dM,Omega=Omega,X=X,Zb=Zb,t=t)

    B    <- bgq$B
    g    <- bgq$g
    q    <- bgq$q
    Xind <- bgq$Xind

    sig   <- Sigma(Z=Z,Zb=Zb,R1=R1,R0=R0,dM=dM,dML=dML,B=B,q=q)
    S     <- sig$S
    psi   <- sig$psi
    eta   <- sig$eta
    n     <- nrow(Z)
    Sbeta <- (n**2)*solve(A)%*%S%*%solve(A)

    m <- length(t)
    y <- cumsum(du[1:m])

#   return value
    obj <- list(beta = beta,
                var  = Sbeta/n,
                y    = y,
                t    = t,
                w    = w,
                varnames = colnames(Z),
                i    = iter,
                call = match.call())
    class(obj) <- "CompoML"


    return(obj)
}



## cox model for censoring

#' @import survival
coxph.LTFU=function(t,Z,X,delta,dNL,Y,YL,ep){

    n=nrow(Z)
    p=ncol(Z)
    m=length(t)

    model=coxph(Surv(X,delta==0)~Z)
    g=summary(model)$coef[,1]
    l=basehaz(model,centered=FALSE)[,1]
    dL.crude=l-c(0,l[-length(l)])


    surv=exp(-basehaz(model,centered=FALSE)[,1]%*%t(exp(Z%*%g)))
    surv=rbind(rep(1,n),surv)


    #fitted.t=c(0,X[order(X)])
    fitted.t=c(0,basehaz(model,centered=FALSE)[,2])
    dL.crude=c(0,dL.crude)
    G=matrix(NA,n,m)
    dL=rep(0,m)

    for (j in 1:m){
        G[,j]=surv[sum(fitted.t<=t[j]),]
        if (any(abs(fitted.t-t[j])<ep)){
            dL[j]=dL.crude[abs(fitted.t-t[j])<ep][1]}
    }
    dML <- dNL-Y*(exp(Z%*%g)%*%t(dL))

    W=matrix(0,n,m)
    for (i in 1:n)
    {
        for (j in 1:m)
        {
            denom <- G[i,sum(min(X[i],t[j])>=t)]
            if (!is.na(denom)& denom!=0)
            {
                W[i,j]=YL[i,j]*G[i,j]/denom
            }
        }
    }
    return(list(W=W,dL=dL,dML=dML,gamma=g))
}

## compute martingale residuals
event.mart=function(beta,W,Z,dN)
{
    p=ncol(Z)
    m=ncol(W)
    n=nrow(Z)
    S0=rep(NA,m)
    S1=list()
    length(S1)=m
    S2=list()
    length(S2)=m
    Zb=matrix(NA,p,m)
    du=rep(NA,m)

    expZbeta <- exp(Z%*%beta)
    for (j in 1:m)
    {
        S0[j]   <- mean(W[,j]*expZbeta)
        S1[[j]] <- t(Z)%*%(W[,j]*expZbeta)/n

        if (S0[j]>0){Zb[,j]=S1[[j]]/S0[j]}else{Zb[,j]=0}
        du[j] <- sum(dN[,j]/S0[j],na.rm=T)/n
        if (is.na(du[j])){du[j]=0}
    }
    dM <- dN - W * (expZbeta %*% t(du))

    return(list(du=du,dM=dM,S0=S0,S1=S1,Zb=Zb))

}

UL <- function(beta,W,Z,dN)
{
    p=ncol(Z)
    m=ncol(W)
    n=nrow(Z)
    U=rep(0,p)
    A=matrix(0,p,p)
    tmp=matrix(0,p,p)
    mart=event.mart(beta,W,Z,dN)
    dM=mart$dM
    du=mart$du
    Zb=mart$Zb
    S0=mart$S0
    S1=mart$S1

    expZbeta <- as.vector(exp(Z %*% beta))

    #S2mat <- Z %*% t(t(expZbeta * W) %*% Z)  / n
    for (j in 1:m)
    {
        U  <- U + (t(Z) - matrix(rep(Zb[,j],n),p,n)) %*% dM[,j]
        #S2 <- (1/n) * t(Z) %*% diag(as.vector(expZbeta * W[,j])) %*% Z
        S2 <- crossprod(Z, as.vector(expZbeta * W[,j]) * Z) / n
        if (S0[j]>0)
        {
            tmp <- tmp + (S2/S0[j]-(S1[[j]]/S0[j]) %*% t(S1[[j]]/S0[j])) * sum(dN[,j])
        }
    }

    A=tmp

    return(list(U=U,A=A))
}


NRbeta=function(W,Z,dN,init,err)
{

    beta <- init
    obj  <- UL(beta = beta,
               W    = W,
               Z    = Z,
               dN   = dN)
    A    <- obj$A
    U    <- obj$U
    inc  <- NA
    if (det(A)!=0)
    {
        inc <- solve(A, U)
    }
    #if (is.na(inc)){inc=0}

    i <- 0
    while(sum(abs(inc)) > err)
    {
        beta <- beta + inc
        obj  <- UL(beta = beta,
                   W    = W,
                   Z    = Z,
                   dN   = dN)
        A    <- obj$A
        U    <- obj$U
        inc  <- solve(A, U)
        i    <- i+1
        #print(U)
        #print(beta)
    }
    return(list(beta = beta, U = U, A = A, i = i))

}



Omega=function(gamma,Y,Z,dNL)
{
    n=nrow(Y)
    p=ncol(Z)
    m=ncol(Y)
    R0=1/n*t(Y)%*%exp(Z%*%gamma)
    R1=1/n*t(Y)%*%(matrix(rep(exp(Z%*%gamma),p),n,p)*Z)

    tmp=matrix(0,p,p)
    expZgamma <- exp(Z%*%gamma)
    for (j in 1:m)
    {

        #R2 <- (1/n) * t(Z) %*% diag(as.vector(expZgamma * Y[,j])) %*% Z
        R2 <- crossprod(Z, as.vector(expZgamma * Y[,j]) * Z) / n
        if (R0[j]>0)
        {
            tmp <- tmp + (R2/R0[j]-(R1[j,]/R0[j])%*%t(R1[j,]/R0[j]))*sum(dNL[,j])
        }
    }
    Omega=tmp/n
    return(list(Omega=Omega,R0=R0,R1=R1))
}


BGQ=function(gamma,Z,R0,R1,dL,dM,Omega,X,Zb,t)
{

    n=nrow(Z)
    p=ncol(Z)
    m=ncol(dM)

    tmp1=matrix(0,m,p)
    tmp2=matrix(0,p,p)
    g=matrix(0,n*m,p)
    Xind=rep(NA,n)
    OmegaInv=matrix(0,p,p)
    if (det(Omega)!=0)
    {
        OmegaInv <- solve(Omega)
    }


    exptgammaZ <- Z %*% exp(gamma)

    for (i in 1:n)
    {
        Xind[i] <- which(t>X[i])[1]
        if (!is.na(Xind[i]))
        {
            exptgammaZi <- exptgammaZ[i] #as.numeric(exp(t(gamma)) %*% Z[i,])
            for (j in Xind[i]:m)
            {
                if (p==1)
                {
                    tmp1[j,] <- tmp1[j,] + exptgammaZi *
                        (t((matrix(rep(Z[i,],m),p,m)-Zb)[,j:m]) %*% dM[i,j:m])
                } else
                {
                    tmp1[j,] <- tmp1[j,] + exptgammaZi *
                        (as.matrix((matrix(rep(Z[i,],m),p,m)-Zb)[,j:m]) %*% dM[i,j:m])
                }
                k  <- j-Xind[i]+1
                Rm <- t(R1[Xind[i]:j,] / matrix(rep(R0[Xind[i]:j],p),k,p))
                Rm[is.na(Rm)] <- 0
                g[(i-1)*m+j,] <- exptgammaZi * as.vector((matrix(rep(Z[i,],k),p,k)-Rm) %*% dL[Xind[i]:j])
                tmp2 <- tmp2 + (Z[i,]-Zb[,j]) %*% t(g[(i-1) * m+j,]) %*% OmegaInv * dM[i,j]
            }
        }
    }
    q <- -tmp1 / n
    B <- -tmp2 / n

    return(list(B=B,q=q,g=g,Xind=Xind))
}


Sigma=function(Z,Zb,R1,R0,dM,dML,B,q)
{
    n=nrow(Z)
    p=ncol(Z)
    m=ncol(dM)

    Rm1 <- t(R1/matrix(rep(R0,p),m,p))
    Rm1[is.na(Rm1)] <- 0
    Rm2 <- t(q/matrix(rep(R0,p),m,p))
    Rm2[is.na(Rm2)] <- 0

    psi=matrix(NA,n,p)
    eta=matrix(NA,n,p)

    for (i in 1:n)
    {
        psi[i,] <- B %*% (matrix(rep(Z[i,],m),p,m)-Rm1) %*% dML[i,] + Rm2 %*% dML[i,]
        eta[i,] <- (matrix(rep(Z[i,],m),p,m)-Zb)%*%dM[i,]
    }

    S <- 1/(n-p-1) * t(psi+eta) %*% (psi+eta)
    #psi=matrix(0,n,p)
    return(list(S=S,psi=psi,eta=eta))
}


Sigmamu=function(Z,Zb,dM,dML,Omega,S0,S1,R0,R1,gamma,du,g,Y,A,eta,psi){
    n=nrow(Z)
    p=ncol(Z)
    m=ncol(dM)

    H=matrix(0,m,p)
    P1=matrix(0,m,m)
    P2=matrix(0,m,p)
    dp1=rep(0,m)
    dp2=matrix(0,m,p)
    phi=matrix(0,n,m)
    invS0=rep(0,m)
    invR0=rep(0,m)
    Ut=matrix(0,p,m)
    dUt=matrix(0,p,m)

    OmegaInv <- solve(Omega)
    for (j in 1:m)
    {
        if (S0[j]>0){invS0[j]=1/S0[j]}
        if (R0[j]>0){invR0[j]=1/R0[j]}
        H[j,]=-Zb[,1:j]%*%du[1:j]
        dUt[,j]=(t(Z)-matrix(rep(Zb[,j],n),p,n))%*%dM[,j]
        #Ut[,j]=rowSums(matrix(as.vector(dUt[,1:j]),p,j))

        for (i in 1:n)
        {
            if (S0[j]>0)
            {
                dp2[j,]=dp2[j,]-1/n*t(g[(i-1)*m+j,]) %*% OmegaInv * (1-Y[i,j])*dM[i,j]/S0[j]
                dp1[j]=dp1[j]-1/n*(1-Y[i,j])*as.numeric(exp(t(gamma)%*%Z[i,]))*dM[i,j]/S0[j]
            }
        }
    }



    Ainv <- solve(A)

    for (j in 1:m)
    {


        P2[j,]=sum(dp2[1:j])
        for (k in 1:j)
        {
            P1[k,j] <- sum(dp1[k:j])
        }
        for (i in 1:n)
        {
            phi1=sum(dM[i,1:j]*invS0[1:j])
            phi2=sum(P1[1:j,j]*dML[i,1:j]*invR0[1:j])
            phi3=t(P2[j,])%*%(matrix(rep(Z[i,],m),p,m)-t(R1*matrix(rep(invR0,p),m,p)))%*%dML[i,]
            #phi4=sqrt(n)*t(H[j,])%*%solve(A)%*%colSums(eta+psi)
            phi4=n*t(H[j,]) %*% Ainv %*% (eta[i,]+psi[i,])
            phi[i,j]=phi1+phi2+phi3+phi4
        }

    }

    xi=colMeans(phi**2)

    return(list(xi=xi,dUt=dUt))
}





