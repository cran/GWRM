\name{GWRM.display}
\Rdversion{1.1}
\alias{GWRM.display}
\title{
Summarizing Generalized Waring Regression Model Fits
}
\description{
Function to produce summaries of the results of \code{GWRM.fit} function.
}
\usage{
GWRM.display(model)
}
\arguments{
  \item{model}{
A model fitted by \code{GWRM.fit} function.
}
}
\value{
\item{Table}{A data frame with the matrix of estimated coefficients of the regressors, standard errors, z-values and p-values.}
\item{betaII}{A data frame with the estimated parameters of the BetaII distribution that models proneness. First, it includes estimates of the parameters that are used in the optimization of the log-likelihood function, \code{betak} and \code{betaro}, such that \code{k=exp(betak)} and \code{ro=1+exp(betaro)}. Then, \code{k} and \code{ro} estimates appear in the third column.}
\item{Fit}{Log-likelihood value, AIC, BIC and degrees of freedom of \code{model}.}
\item{Convergence}{A code that indicates successful convergence of the fitter function used (see \code{nlm} or \code{optim} helps).}
\item{Method}{The code of the fitter function used. See \code{GWRM.fit} help.}
}
\references{
 Rodriguez-Avi, J; Conde-Sanchez, A; Saez-Castillo, A.J., Olmo-Jimenez, M. J. and Martinez Rodriguez, A. M.(2009). A generalized Waring regression model for count data. Computational Statistics and Data Analysis, 53, pp. 3717-3725.
}
\author{
A.J. Saez-Castillo \email{ajsaez@ujaen.es}, J. Rodriguez-Avi, A. Conde-Sanchez, M.J. Olmo-Jimenez, A.M. Martinez-Rodriguez.
}
\examples{
data(goals)

fit<-GWRM.fit(goals~clasif+position+log(played)+offset(played),data=goals)
GWRM.display(fit)
}
\keyword{regression}
\keyword{models}