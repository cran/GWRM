\name{add1}
\alias{add1.gw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add All Possible Single Terms to a GWRM Model
}
\description{
Compute all the single terms in the scope argument that can be added to the GWRM model, fit those models and compute a table of the changes in fit.
}
\usage{
\method{add1}{gw}(object, scope, test = c("none", "Chisq"), k = 2, trace = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted \code{"gw"} object.
}
  \item{scope}{
a formula giving the terms to be considered for adding.
}
\item{test}{
\code{"none"}, which considers the AIC criterion, or \code{Chisq}, which is the likelihood-ratio test.
}
  \item{k}{
the penalty constant in AIC / Cp.
}
  \item{trace}{
if \code{TRUE}, print out progress reports.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
The method is commonly used by \code{\link{step}}.

In a scope formula \code{.} means 'what is already there'.

The default output table gives AIC, defined as minus twice log likelihood plus \code{2p} where \code{p} is the rank of the model (the number of effective parameters).
}
\value{
An object of class "anova" summarizing the differences in fit between the models.
}
\author{
A.J. Saez-Castillo <ajsaez@ujaen.es>, S. Vilchez-Lopez, J. Rodriguez-Avi, A. Conde-Sanchez, M.J. Olmo-Jimenez, A.M. Martinez-Rodriguez

The design was inspired by the R function of the same name.
}
\seealso{
\code{\link{step}}, \code{\link{drop1.gw}}, \code{\link{gw}}
}

\examples{
data(goals)

fit0 <- gw(goals ~ offset(log(played)), data = goals)
summary(fit0)
           
fit1 <- step(fit0, ~ position)
summary(fit1)
}