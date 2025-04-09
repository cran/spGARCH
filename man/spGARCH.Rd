\name{spGARCH}
\alias{spGARCH}
\title{spGARCH - Package to fit spatial ARCH models}

\description{
A collection of functions for simulating and fitting spatial autoregressive conditional heteroscedasticity (spARCH) processes are provided.

The functions \code{\link[spGARCH]{sim.spARCH}} and \code{\link[spGARCH]{sim.spGARCH}} are the main function for simulating spARCH and spGARCH processes, respectively. Via the argument \code{type} several types of spatial ARCH and GARCH can be simulated, e.g., exponential spARCH models, spARCH for oriented processes, or spARCH processes with truncated error support. For details, refer the paper Otto, Schmid, and Garthoff (2018), see \doi{10.1016/j.spasta.2018.07.005}.
Moreover, the package provides function for fitting spARCH models. Basically, there are two functions to fit these kind of model: \code{\link[spGARCH]{qml.spARCH}} and \code{\link[spGARCH]{qml.SARspARCH}}. First, spARCH models can be fitted by \code{\link[spGARCH]{qml.spARCH}}.
}
\author{
Philipp Otto \email{philipp.otto@glasgow.ac.uk}
}
\references{
Philipp Otto, Wolfgang Schmid, Robert Garthoff (2018). Generalised Spatial and Spatiotemporal Autoregressive Conditional Heteroscedasticity. Spatial Statistics 26, pp. 125-145. \doi{10.1016/j.spasta.2018.07.005}, arXiv: \doi{10.48550/arXiv.1609.00711}
}
