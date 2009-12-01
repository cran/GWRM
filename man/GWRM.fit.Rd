\name{GWRM.fit}
\Rdversion{1.1}
\alias{GWRM.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Generalized Waring Regression Models
}
\description{
Function to fit a Generalized Waring Regression Model (GWRM), which is specified giving the response variable and a symbolic description of the linear predictor as a formula.
}
\usage{
GWRM.fit(formula, f = NULL, p0k = 0, p0ro = 1, p0beta = NULL, 
iters = 10000, data, method = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
An object of class "formula": a symbolic description of the model to be fitted.
}
  \item{f}{
An optional vector of frequencies or prior weights to be used in the fitting process. It should be NULL or a numeric vector.
}
  \item{p0k}{
A starting value for the parameter \code{k} in the Generalized Waring distribution.
}
  \item{p0ro}{
A starting value for the parameter \code{ro} in the Generalized Waring distribution.
}
  \item{p0beta}{
A vector of starting values for the coefficients of the regressors.
}
  \item{iters}{
A positive integer specifying the maximum number of iterations to be performed in the optimization process of the log-likelihood.
}
  \item{data}{
A data frame containing the variables in \code{formula}.
}
  \item{method}{
An integer specifying the method to be used for the log-likelihood optimization. See "Details".
}
}
\details{
Method 1 indicates optimization with \code{nlm} function. Method 2 to 6 indicate optimization with "Nelder-Mead", "BFGS", "CG", "L-BFGS-B" and "SANN" methods of \code{optim} function respectively. 

An offset may be included in \code{formula}. GWRM is a multiplicative model, so the offset is added as a multiplicative term. For example, if \code{formula=y~x+offset(z)} then, the model is \code{y=z*exp(beta0+beta1*x)}.
}
\value{
\item{dataset}{A data frame with data and vector \code{f} of frequencies or weights.}
\item{response}{The response variable.}
\item{model}{The formula for the model.}
\item{covars}{The names of the predictors.}
\item{offset}{A logical value indicating whether \code{formula} includes an offset.}
\item{optimum}{The optimum of log-likelihood function.}
\item{aic}{Akaike Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters.}
\item{bic}{Bayesian Information Criterion, minus twice the maximized log-likelihood plus the number of parameters multiplied by the logarithm of the number of observations.}
\item{df}{The degrees of freedom of the model.}
\item{coefficients}{The estimates of the coefficients of the model. The last two correspond to the parameters of the BetaII distribution. It should be taken into account that they are parametrized as \code{k=exp(betak)} and \code{ro=1+exp(betaro)}; here you'll find \code{betak} and \code{betaro}.}
\item{betaIIpars}{The estimates of the parameters \code{k} and \code{ro} of the BetaII distribution that models proneness.}
\item{betascoefs}{The estimates of the coefficients of the linear predictor.}
\item{fitted.values}{The fitted mean values.}
\item{hessian}{A symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function.}
\item{cov}{A matrix giving an estimate of the covariance matrix of the coefficients of the model.}
\item{se}{A vector of estimates of the standard errors of the estimated coefficients.}
\item{corr}{A matrix giving an estimate of the correlation matrix of the coefficients of the model.}
\item{code}{A code that indicates successful convergence of the fitter function used (see \code{nlm} and \code{optim} helps).}
\item{method}{The name of the fitter function used.}
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
