\name{partvar}
\alias{partvar}
\alias{partvar.default}
\title{
Partition of the variance in a GWRM model
}
\description{
In a GWRM model, the variance may be split into three terms. The first component of this decomposition represents the variability due to randomness and it comes from the underlying Poisson model. The other two components refer to the variability that is not due to randomness but is explained by the presence of liability and proneness, respectively.
}
\usage{
\method{partvar}{default}(object, newdata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object class \code{"gw"} for which the partition is desired.
}
  \item{newdata}{
optionally, a data frame in which to look for variables with which to obtain the partition. If omitted, all the cases are used.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
One of the main drawbacks of using the Univariate Generalized Waring Distribution with parameters \emph{a}, \emph{k} and \emph{ro} is that the parameters \emph{a} and \emph{k} are interchangeable when there is no auxiliary information given by the covariates. This identification problem prevents liability and proneness components from being distinguished in the univariate fits. To solve it, Irwin (1968) proposed that the expert should deduce which of these components is which from their own knowledge of the phenomenon. Xekalaki (1984) proposed a less subjective solution, developing a bivariate model that divides the observation period into two non-overlapping subperiods in which the model for proneness does not change. In the GWRM with, at least, one covariate, the parameters \emph{a} and \emph{k} are not interchangeable because, as in Xekalaki's bivariate model, the random model for proneness does not change. So, the identification problem of the non-random components is solved.
}
\value{
Two data frames, with ratio of sources of variation and sources of variation in which variance is splitted.
}
\references{
J.O. Irwin (1968). The generalized waring distribution applied to accident theory Journal of the Royal Statistical Society. Series A, 131 (2), pp. 205-225

Rodriguez-Avi, J., Conde-Sanchez, A., Saez-Castillo, A. J., Olmo-Jimenez, M. J. and Martinez Rodriguez, A. M.(2009). A generalized Waring regression model for count data. Computational Statistics and Data Analysis, 53, pp. 3717-3725

E. Xekalaki (1984). The bivariate generalized waring distribution and its application to accident theory Journal of the Royal Statistical Society. Series A, 147 (3), pp. 488-498
}
\author{
A.J. Saez-Castillo <ajsaez@ujaen.es>, S. Vilchez-Lopez, J. Rodriguez-Avi, A. Conde-Sanchez, M.J. Olmo-Jimenez, A.M. Martinez-Rodriguez
}

\seealso{
\code{\link{gw}}, \code{\link{predict.gw}}
}

\examples{
data(goals)

fit <- gw(goals ~ position + offset(log(played)), data = goals)

pos <- factor(c("Defender", "Midfielder"), levels = c("Defender", "Midfielder", "Forward"))
lev <- data.frame(position = pos, played = c(17, 21))

partvar(fit, newdata = lev)
}

