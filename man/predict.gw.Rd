\name{predict}
\alias{predict.gw}
\title{
GWRM Predictions
}
\description{
Obtains predictions from a fitted GWRM object.
}
\usage{
\method{predict}{gw}(object = NULL, newdata = NULL, ...)
}
\arguments{
  \item{object}{
a fitted object of class inheriting from \code{"gw"}.
}
  \item{newdata}{
optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\value{
A data frame with \code{newdata} and their fitted means.
}
\seealso{
\code{\link{gw}}, \code{\link{partvar}}
}
\examples{
data(goals)

fit <- gw(goals ~ position + offset(log(played)), data = goals)

pos <- factor(c("Defender", "Midfielder"), levels = c("Defender", "Midfielder", "Forward"))
lev <- data.frame(position = pos, played = c(17, 21))

predict(fit, newdata = lev)
}