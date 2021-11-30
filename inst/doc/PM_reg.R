## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=F-------------------------------------------------------------------
#  obj <- CompoML(id,time,status,Z)

## ----eval=F-------------------------------------------------------------------
#  plot(obj, z)

## ----setup--------------------------------------------------------------------
## load the package and data
library(Wcompo)
data(hfmock)
head(hfmock)

## -----------------------------------------------------------------------------
## fit a weighted PM (w_D=2, w_1=1)
obj <- CompoML(hfmock$id,hfmock$time,hfmock$status,hfmock[,c("Training","HF.etiology")],
               w=c(2,1))
## print out the result
obj

## ---- fig.height = 4.5, fig.width=7.2-----------------------------------------
oldpar <- par(mfrow = par("mfrow"))
par(mfrow=c(1,2))
## plot the estimated mean function for
## non-ischemic patients by treatment
plot(obj,c(1,0),ylim=c(0,1.5),xlim=c(0,50),
     main="Non-ischemic",
     xlab="Time (months)",cex.main=1.2,lwd=2)
plot(obj,c(0,0),add=TRUE,cex.main=1.2,lwd=2,lty=2)
legend("topleft",lty=1:2,lwd=2,c("Exercise training","Usual care"))


## plot the estimated mean function for
## ischemic patients by treatment
plot(obj,c(1,1),ylim=c(0,1.5),xlim=c(0,50),
     main="Ischemic",
     xlab="Time (months)",cex.main=1.2,lwd=2)
plot(obj,c(0,1),add=TRUE,cex.main=1.2,lwd=2,lty=2)
legend("topleft",lty=1:2,lwd=2,c("Exercise training","Usual care"))
par(oldpar)

