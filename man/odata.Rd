\name{odata}
\alias{odata}
\docType{data}
\title{
Case-control dataset
}
\description{
An example of the case-control data (or the outcome data).
}
\usage{data("odata")}
\format{
  A data frame with 3500 observations on the following 7 variables.
  \describe{
    \item{\code{d}}{case-control status}
    \item{\code{x1}}{numeric covariate}
    \item{\code{x2}}{character covariate}
    \item{\code{g1}}{instrument}
    \item{\code{g2}}{instrument}
    \item{\code{g3}}{instrument}
    \item{\code{id}}{subject IDs}
  }
}
\details{
The column \code{id} is always required in exposure data.
}

\examples{
data(odata)
head(odata, 2)
}
\keyword{datasets}
