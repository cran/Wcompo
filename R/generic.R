
##generic functions


#' Print the analysis results of the proportional means model
#'
#' Print the analysis results of the proportional means model.
#'
#' @param x An object returned by \code{\link{CompoML}}.
#' @param ... Further arguments passed to or from other methods.
#' @return Print the results of \code{\link{CompoML}} object
#' @importFrom stats pnorm printCoefmat qnorm
#' @keywords CompoML
#' @seealso \code{\link{CompoML}}, \code{\link{plot.CompoML}}.
#' @export
print.CompoML=function(x,...){
cat("\n")
cat("Call:\n")
print(x$call)
cat("\n")

  beta=x$beta
  var=x$var
  t=x$t
  y=x$y


  w=x$w

  varnames=x$varnames

  se=sqrt(diag(var))

  pv=2*(1-pnorm(abs(beta/se)))


cat("Proportional means regression models for composite endpoint\n (Mao and Lin, 2016, Biostatistics):\n\n")

K=length(w)
weight=matrix(w,1,K)
rownames(weight)="Weight"
colnames(weight)=rep(" ",K)
colnames(weight)[1]="Event 1 (Death)"
for (k in 2:K){
colnames(weight)[k]=paste("Event",k)
}
print(weight)
cat("\n\n")
cat("Newton-Raphson algorithm converges in", x$i,"iterations.\n\n")


table=cbind(Estimate = beta,
             StdErr = se,
             z.value = beta/se,
             p.value = pv)

colnames(table)=c("Estimate","se","z.value","p.value")
rownames(table)=varnames

cat("Estimates for Regression parameters:\n\n")
printCoefmat(table, P.values=TRUE, has.Pvalue=TRUE)
cat("\n")
cat("\n")

za=qnorm(0.975)
MR=cbind(exp(beta),exp(beta-za*se),exp(beta+za*se))
colnames(MR)=c("Mean Ratio","95% lower CL","95% higher CL")
rownames(MR)=varnames
cat("Point and interval estimates for the mean ratios:\n\n")
print(MR)
cat("\n")

}




#' Plot the predicted mean function under the proportional means model
#'
#' Plot the predicted mean function under the proportional means model for
#' a new observation.
#'
#' @param x An object returned by \code{\link{CompoML}}.
#' @param z Covariate vector for the new observation. If \code{NULL}, the baseline
#' mean function will be plotted.
#' @param xlab A label for the x axis.
#' @param ylab A label for the y axis.
#' @param lty Line type for the plot.
#' @param frame.plot Boolean argument indicating whether to add a rectangular frame
#' to the plot.
#' @param add If TRUE, the curve will be overlaid on an existing plot; otherwise,
#' a separate plot will be constructed.
#' @param ... Other arguments that can be passed to the underlying \code{plot} method.
#' @return No return value, called for side effects.
#' @seealso \code{\link{CompoML}}, \code{\link{print.CompoML}}.
#' @keywords CompoML
#' @importFrom graphics lines
#' @export
#' @examples
#' ## see example for CompoML
plot.CompoML=function(x,z=NULL,xlab="Time",
ylab="Mean function",lty=1,frame.plot=FALSE,add=FALSE,
...){

  beta=x$beta
  var=x$var
  t=c(0,x$t)
  y=x$y

  p<-length(beta)

  if(is.null(z)){
    z<-rep()
  }


  mf=c(0,exp(sum(beta*z))*y)



  if (add==F){
	plot(t,mf, type="l", lty=lty,xlab=xlab,ylab=ylab,frame.plot=frame.plot,...)
  }else{
    lines(t,mf,lty=lty,...)
  }
  }

