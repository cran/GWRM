\name{GWRM.stats}
\Rdversion{1.1}
\alias{GWRM.stats}
\title{
Summary of the Generalized Waring Distribution of the cases specified
}
\description{
In a Generalized Waring Regression Model (GWRM) with parameters \code{k} and \code{ro}, the probability distribution of each case is a Univariate Generalized Waring Distribution (UGWD) with parameters \code{a}, \code{k} and \code{ro}, where \code{a} is determined by the covariates values. This function calculates, for the cases specified by the arguments \code{covars} in \code{model}, the value of \code{a}, the mean, a confidence interval for the mean, the variance and the partition of the variance in terms of randomness, liability and proneness.
}
\usage{
GWRM.stats(model = NULL, covars = NULL, alpha = 0.05)
}
\arguments{
  \item{model}{
A model given by \code{GWRM.fit} function.
}
  \item{covars}{
A data frame whose columns must have the same names than the covariates of the GWRM specified by the argument \code{model}. Rows specify each case that you wish to analyze.
}
  \item{alpha}{
The confidence level intervals for the means.
}
}
\details{
If any of the covariates is a factor, it must be introduced in \code{covars} also as a factor.
}
\value{
\item{params}{A data frame with the estimated parameters corresponding to the UGWD of each case specified by \code{covars} in the GWRM indicated in \code{model}.}
\item{stats}{A data frame with the estimated mean, the confidence interval of the mean and the variance of the UGWD of each case specified by \code{covars} in the GWRM indicated in \code{model}.}
\item{partvar}{A data frame with the partition of the variance (in rates) of the UGWD of each case specified by \code{covars} in the GWRM indicated in \code{model}; the three components of the variance are randomness, liability and proneness.}
}
\references{
Rodriguez-Avi, J; Conde-Sanchez, A; Saez-Castillo, A.J., Olmo-Jimenez, M. J. and Martinez Rodriguez, A. M.(2009). A generalized Waring regression model for count data. Computational Statistics and Data Analysis, 53, pp. 3717-3725.
}
\author{
A.J. Saez-Castillo \email{ajsaez@ujaen.es}, J. Rodriguez-Avi, A. Conde-Sanchez, M.J. Olmo-Jimenez, A.M. Martinez-Rodriguez
}
\examples{
data(goals)

fit<-GWRM.fit(goals~clasif+position+log(played)+offset(played),data=goals)
GWRM.display(fit)

dem<-factor(c("Forward","Defender"),levels=c("Defender","Midfielder","Forward"))
niveles<-data.frame(clasif=c(1,4),played=c(12,15),position=dem)
GWRM.stats(model=fit,covars=niveles)
}
\keyword{regression}
\keyword{models}