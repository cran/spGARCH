\name{extractAIC.spARCH}
\alias{extractAIC.spARCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the Akaike information criterion
}
\description{
The function extracts the log-likelihood of a spatial ARCH model.
}
\usage{
## S3 method for class 'spARCH'
\method{extractAIC}{spARCH}(fit, scale, k = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
\code{spARCH} object (fitted model resulting of \code{qml.spARCH}).
}
  \item{scale}{
currently unused for \code{spARCH} objects.
}
\item{k}{
parameter specifying the weight for the penalizing term.
}
  \item{...}{
Other arguments.
}
}
\details{
Numeric vector of length 2 is returned. The first element specify the \code{edf} (equivalent degree of freedom) and the Akaike information criterion is returned as second element.
}
\author{
Philipp Otto \email{potto@europa-uni.de}
}
