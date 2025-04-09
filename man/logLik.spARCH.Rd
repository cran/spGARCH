\name{logLik.spARCH}
\alias{logLik.spARCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract logarithmic likelihood
}
\description{
The function extracts the log-likelihood of a spatial ARCH model.
}
\usage{
## S3 method for class 'spARCH'
\method{logLik}{spARCH}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
\code{spARCH} Object of \code{qml.spARCH}.
}
  \item{...}{
Other arguments.
}
}
\details{
\code{logLik} object is returned.
}
\author{
Philipp Otto \email{philipp.otto@glasgow.ac.uk}
}
\seealso{
\code{\link[spGARCH]{residuals.spARCH}}, \code{\link[spGARCH]{fitted.spARCH}}.
}
