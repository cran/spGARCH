\name{spGARCH}
\alias{spGARCH}
\title{spGARCH - Package to fit spatial ARCH models}

\description{
A collection of functions for simulating and fitting spatial autoregressive conditional heteroscedasticity (spARCH) processes are provided.

The function \code{\link{sim.spARCH}} is the main function for simulating spARCH processes. Via the argument \code{type} several types of spatial ARCH can be simulated, e.g., exponential spARCH models, spARCH for oriented processes, or spARCH processes with truncated error support. For details, refer the paper Otto, Schmid, and Garthoff (2017).

Moreover, the package provides function for fitting spARCH models. Basically, there are two functions to fit these kind of model: \code{\link{qml.spARCH}} and \code{\link{qml.SARspARCH}}. First, spARCH models can be fitted by \code{\link{qml.spARCH}}.
}
\author{
Philipp Otto \email{potto@europa-uni.de}
}
\references{
Philipp Otto, Wolfgang Schmid, Robert Garthoff (20xx). Spatial and Spatiotemporal Autoregressive Conditional Heteroscedasticity. \emph{Journal ...}
}
