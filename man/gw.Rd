\name{gw}
\alias{gw}
\title{
Fitting GWRM Models
}
\description{
\code{gw} is used to fit Generalized Waring Regression Models (GWRM), specified by giving a symbolic description of the linear predictor.
}
\usage{
gw(formula, data, weights, k=NULL, subset, na.action, kstart=1, rostart=2,
    betastart=NULL, offset, control = list(...), method = NULL, hessian=TRUE,
    model = TRUE, x = FALSE, y = TRUE, ...)
}
\arguments{
  \item{formula}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).
}
  \item{weights}{
an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.
}
  \item{k}{
optional value for the \emph{k} parameter. If \code{NULL}, it is estimated.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{
a function which indicates what should happen when the data contain \code{NAs}. See \code{\link{glm}}.
}
  \item{kstart}{
starting value for the \emph{k} parameter.
}
  \item{rostart}{
starting value for the \emph{ro} parameter.
}
  \item{betastart}{
starting values for the vector of means.
}
  \item{offset}{
this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
}
  \item{control}{
a list of parameters for controlling the fitting process.
}
  \item{method}{
the method to be used in fitting the model. The default method initially uses non-linear minimization (\code{nlm}) and Nelder-Mead optimization (\code{optim}) to fit a model which is then re-fitted by "L-BFGS-B" (\code{optim}). In this way, SE estimates for all the model parameters are provided. \code{"nlm"} and \code{"Nelder-Mead"} are also possible values, but they do not provide SE estimates for \code{k} and \code{ro}.
}
  \item{hessian}{
if \code{TRUE}, the hessian of \emph{f} at the minimum is returned.
}
  \item{model}{
a logical value indicating whether \emph{model frame} should be included as a component of the returned value.
}
  \item{x, y}{
logical values indicating whether the response vector and model matrix used in the fitting process should be returned as components of the returned value.
}
  \item{\dots}{
further arguments.
}
}
\details{
The GWRM is fully described in Rodriguez-Avi et al. (2009).
}
\value{
\code{gw} returns an object of class \code{"gw"}. The function \code{summary} can be used to obtain or print a summary of the results. An object of class \code{"gw"} is a list containing the following components:

\item{Y}{if requested (the default), the \emph{y} vector used.}
\item{W}{the weights supplied, a vector of \code{1}s if none were.}
\item{covars}{names of the covariates in the model.}
\item{nobs}{number of observations.}
\item{covoffset}{a logical value specifying if an offset is present.}
\item{loglik}{the maximized log-likelihood.}
\item{aic}{a version of Akaike's \emph{An Information Criterion}, minus twice the maximized log-likelihood plus twice the number of parameters.}
\item{bic}{Bayesian Information Criterion, minus twice the maximized log-likelihood plus the number of parameters multiplied by the logarithm of the number of observations.}
\item{df.residual}{the residual degrees of freedom.}
\item{residuals}{the residuals in the final iteration of the fit.}
\item{coefficients}{a named vector of coefficients.}
\item{betaIIpars}{parameters estimates of the BetaII distribution.}
\item{betascoefs}{a vector of coefficients.}
\item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
\item{hessian}{a symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function.}
\item{cov}{an estimate of the covariance matrix of the model coefficients.
}
\item{se}{a vector of the standard errors estimates of the estimated coefficients.}
\item{corr}{an estimate of the correlation matrix of the model coefficients.}
\item{code}{a code that indicates successful convergence of the fitter function used (see nlm
and optim helps).}
\item{method}{the name of the fitter function used.}
\item{k}{if requested, the \emph{k} value used.}
\item{kBool}{a logical value specifying whether there is a \emph{k} value or it is estimated.}
\item{call}{the matched call.}
\item{formula}{the formula supplied.}
\item{terms}{the \code{\link{terms}} object used.}
\item{data}{the \code{data} argument.}
\item{offset}{the offset vector used.}
\item{control}{the value of the \code{control} argument used.}
\item{method}{the name of the fitter function used.}
\item{contrasts}{(where relevant) the contrasts used.}
\item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
}
\references{
Rodriguez-Avi, J; Conde-Sanchez, A; Saez-Castillo, A.J., Olmo-Jimenez, M. J. and Martinez Rodriguez, A. M.(2009). A generalized Waring regression model for count data. Computational Statistics and Data Analysis, 53, pp. 3717-3725.
}
\author{
A.J. Saez-Castillo <ajsaez@ujaen.es>, S. Vilchez-Lopez, J. Rodriguez-Avi, A. Conde-Sanchez, M.J. Olmo-Jimenez, A.M. Martinez-Rodriguez

The design was inspired by the \code{\link{glm}} R function.
}

\seealso{
\code{\link{residuals.gw}}, \code{\link{partvar}} 
}
\examples{
data(goals)

fit <- gw(goals ~ position + offset(log(played)), data = goals)
summary(fit)

}

