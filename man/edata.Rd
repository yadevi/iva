\name{edata}
\alias{edata}
\docType{data}
\title{
Exposure dataset
}
\description{
An example of the exposure data.
}
\usage{data("edata")}
\format{
  A data frame with 1000 observations on the following 8 variables.
  \describe{
    \item{\code{d}}{condition status (e.g. disease or not). }
    \item{\code{z}}{quantitative exposure. }
    \item{\code{x2}}{character covariate}
    \item{\code{x3}}{binary covariate}
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
data(edata)
head(edata, 2)
}
\keyword{datasets}

