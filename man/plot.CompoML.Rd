% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R
\name{plot.CompoML}
\alias{plot.CompoML}
\title{Plot the predicted mean function under the proportional means model}
\usage{
\method{plot}{CompoML}(
  x,
  z = NULL,
  xlab = "Time",
  ylab = "Mean function",
  lty = 1,
  frame.plot = FALSE,
  add = FALSE,
  ...
)
}
\arguments{
\item{x}{An object returned by \code{\link{CompoML}}.}

\item{z}{Covariate vector for the new observation. If \code{NULL}, the baseline
mean function will be plotted.}

\item{xlab}{A label for the x axis.}

\item{ylab}{A label for the y axis.}

\item{lty}{Line type for the plot.}

\item{frame.plot}{Boolean argument indicating whether to add a rectangular frame
to the plot.}

\item{add}{If TRUE, the curve will be overlaid on an existing plot; otherwise,
a separate plot will be constructed.}

\item{...}{Other arguments that can be passed to the underlying \code{plot} method.}
}
\value{
No return value, called for side effects.
}
\description{
Plot the predicted mean function under the proportional means model for
a new observation.
}
\examples{
## see example for CompoML
}
\seealso{
\code{\link{CompoML}}, \code{\link{print.CompoML}}.
}
\keyword{CompoML}
