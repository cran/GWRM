\name{GWRM.add}
\Rdversion{1.1}
\alias{GWRM.add}
\title{
Add All Possible Single Terms to a Generalized Waring Regression Model (GWRM)
}
\description{
Computes all the single terms in the scope argument that can be added to the model, fits those models and calculates a table of the changes in fit. 
}
\usage{
GWRM.add(object = NULL, scope = NULL, method=1, iters = 10000)
}
\arguments{
  \item{object}{
A fitted model object obtained with \code{GWRM.fit}.
}
  \item{scope}{
A formula giving the terms to be considered for adding.
}
  \item{iters}{
A positive integer specifying the maximum number of iterations to be performed in the optimization process of the log-likelihood function.
}
  \item{method}{
An integer specifying the method to be used for the log-likelihood optimization. 1 indicates optimization with \code{nlm} function; 2 to 6 indicate optimization with "Nelder-Mead", "BFGS", "CG", "L-BFGS-B" and "SANN" methods of \code{optim} function respectively.   
}
}
\value{
A data frame with the names of the covariates which are added to the model, the value of the chi-squared statistic in the comparison of the new models with \code{object} argument, the degrees of freedom and the associated p-value.
}
\references{
 Rodriguez-Avi, J; Conde-Sanchez, A; Saez-Castillo, A.J., Olmo-Jimenez, M. J. and Martinez Rodriguez, A. M.(2009). A generalized Waring regression model for count data. Computational Statistics and Data Analysis, 53, pp. 3717-3725.
}
\author{
A.J. Saez-Castillo \email{ajsaez@ujaen.es}, J. Rodriguez-Avi, A. Conde-Sanchez, M.J. Olmo-Jimenez, A.M. Martinez-Rodriguez
}
\examples{
data(goals)

fit.0<-GWRM.fit(goals~offset(played),data=goals)
GWRM.display(fit.0)
GWRM.add(object=fit.0,scope=~clasif+position)
}
\keyword{regression}
\keyword{models}