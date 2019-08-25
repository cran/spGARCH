\name{spGARCH}
\alias{spGARCH}
\title{spGARCH - Package to fit spatial ARCH models}

\description{
A collection of functions for simulating and fitting spatial autoregressive conditional heteroscedasticity (spARCH) processes are provided.

The functions \code{\link{sim.spARCH}} and \code{\link{sim.spGARCH}} are the main function for simulating spARCH and spGARCH processes, respectively. Via the argument \code{type} several types of spatial ARCH and GARCH can be simulated, e.g., exponential spARCH models, spARCH for oriented processes, or spARCH processes with truncated error support. For details, refer the paper Otto, Schmid, and Garthoff (2018) published in \href{https://doi.org/10.1016/j.spasta.2018.07.005}{Spatial Statistics}.

Moreover, the package provides function for fitting spARCH models. Basically, there are two functions to fit these kind of model: \code{\link{qml.spARCH}} and \code{\link{qml.SARspARCH}}. First, spARCH models can be fitted by \code{\link{qml.spARCH}}.
}
\author{
Philipp Otto \email{potto@europa-uni.de}
}
\references{
Philipp Otto, Wolfgang Schmid, Robert Garthoff (2018). Generalised Spatial and Spatiotemporal Autoregressive Conditional Heteroscedasticity. \href{https://doi.org/10.1016/j.spasta.2018.07.005}{Spatial Statistics} 26, pp. 125-145.
}
