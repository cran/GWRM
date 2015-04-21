\name{residuals}
\alias{residuals.gw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract and Visualize Model Residuals}

\description{
\code{residuals} is a method which extracts model residuals from \code{"gw"}, commonly returned by \code{gw} function. Optionally, it produces a normal plot with simulated envelope of the residual.
}
\usage{
\method{residuals}{gw}(object, type = "pearson", rep = 19, 
      envelope = TRUE, title = "Simulated Envelope of Residuals",trace = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class \code{"gw"} holding the fitted model
}
  \item{type}{
type of residuals to be extracted. Default is \code{pearson}. \code{response} and \code{deviance} are also available.
}
  \item{rep}{
number of replications for envelope construction. Default is 19, that is the smallest 95\% band that can be build
}
  \item{envelope}{
a logical value to specify if the envelope is required.
}
  \item{title}{
a title for the envelope.
}
  \item{trace}{
if \code{TRUE} a sort of information is printed during the running time.
}
  \item{\dots}{
further arguments.
}
}
\details{
The usual Q-Q plot may show an unsatisfactory pattern of the residuals of a model fitted: then we are led to think that the model is badly specificated. The normal plot with simulated envelope indicates that under the distribution of the response variable the model is OK if only a few points fall off the envelope.
}
\value{
\item{type}{residuals type (\code{"pearson"}, \code{"deviance"} or \code{"response"}) required}
\item{residuals}{actual residuals of the model}
\item{sim.residuals}{simulated residuals}
}

\author{
A.J. Saez-Castillo <ajsaez@ujaen.es>, S. Vilchez-Lopez, J. Rodriguez-Avi, A. Conde-Sanchez, M.J. Olmo-Jimenez, A.M. Martinez-Rodriguez

The design was inspired by the \code{envelope.pgam} R function of \code{pgam} package.
}

\examples{
data(goals)
fit <- gw(goals ~ position + offset(log(played)), data = goals)
res <- residuals(fit, envelope = FALSE, type = "deviance")
normal.score <- qnorm(1:fit$nobs / (fit$nobs + 1))
plot(normal.score, res$residuals, type = "l", 
  xlab = "Standard normal quantiles", 
  ylab = "Residuals (deviance)")
}