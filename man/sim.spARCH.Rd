\name{sim.spARCH}
\alias{sim.spARCH}
\title{Simulation of spatial ARCH models}
\description{
The function generates \code{n} random numbers of a spatial ARCH process for given parameters and weighting schemes.
}
\usage{
sim.spARCH(n = dim(W)[1], rho, alpha, W, b = 2, type = "gaussian", control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of observations. If \code{length(n)} \eqn{> 1}{> 1}, the length is taken to be the number required. Default \code{dim(W)[1]}}
  \item{rho}{spatial dependence parameter rho}
  \item{alpha}{unconditional variance level alpha}
  \item{W}{\code{n} times \code{n} spatial weight matrix}
  \item{b}{parameter \code{b} for spatial exponential ARCH (only needed if \code{type = "exp"}). Default 2.}
  \item{type}{type of simulated spARCH process (see details)}
  \item{control}{list of control arguments (see below)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function simulates \code{n} observations \eqn{Y = (Y_1, ..., Y_n)'}{Y = (Y_1, ..., Y_n)'} of a spatial ARCH process, i.e.,
\deqn{
\boldsymbol{Y} = diag(\boldsymbol{h})^{1/2} \boldsymbol{\varepsilon} \, ,
}{
Y = diag(h)^(1/2) \epsilon ,
}
where \eqn{\boldsymbol{\varepsilon}}{\epsilon} is a spatial White Noise process. The definition of \eqn{\boldsymbol{h}}{h} depends on the chosen \code{type}. The following types are available.
\itemize{
\item \code{type = "gaussian"} - simulates \eqn{\boldsymbol{\varepsilon}}{\epsilon} from a truncated normal distribution on the interval \eqn{[-a, a]}{[-a, a]}, such that \eqn{\boldsymbol{h} > 0}{h > 0} with
\deqn{
\boldsymbol{h} =  \alpha + \rho \mathbf{W} \boldsymbol{Y}^{(2)} \; \mbox{and} \; a = 1 / ||\rho^2\mathbf{W}^2||_1^{1/4}.
}{
h = \alpha + \rho W Y^(2)  and  a = 1 / (\rho^2||W^2||_1)^(1/4)
}
Note that the normal distribution is not trunctated (\eqn{a = \infty}{a =  \infty}), if \eqn{\mathbf{W}}{W} is a strictly triangular matrix, as it is ensured that \eqn{\boldsymbol{h} > \boldsymbol{0}}{h > 0}. Generally, it is sufficient that if there exists a permutation such that \eqn{\mathbf{W}}{W} is strictly triangular. In this case, the process is called oriented spARCH process.
\item \code{type = "exp"} - simulates an exponential spARCH process (E-spARCH), i.e.,
\deqn{
\ln\boldsymbol{h} = \alpha + \rho \mathbf{W} g(\boldsymbol{\varepsilon}) \, .
}{
ln h = \alpha + \rho W g_b(\epsilon) .
}
For the E-spARCH process, the errors follow a standard normal distribution. The function \eqn{g_b}{g_b} is given by
\deqn{
g_b(\boldsymbol{\varepsilon}) = (\ln|\varepsilon(\boldsymbol{s}_1)|^{b}, \ldots, \ln|\varepsilon(\boldsymbol{s}_n)|^{b})' \, .
}{
g_b(\epsilon) = (ln|h(\epsilon_1)|^b, ..., ln|h(\epsilon_n)|^b)' .
}
\item \code{type = "complex"} - allows for complex solutions of \eqn{\boldsymbol{h}^{1/2}}{h^(1/2)} with
\deqn{
\boldsymbol{h} =  \alpha + \rho \mathbf{W} \boldsymbol{Y}^{(2)}  \, .
}{
h = | \alpha + \rho W Y^(2) |.
}
The errors follow a standard normal distribution.
}
}
\value{
The functions returns a vector \eqn{\boldsymbol{y}}{y}.
}
\references{
Philipp Otto, Wolfgang Schmid, Robert Garthoff (2017). Spatial and Spatiotemporal Autoregressive Conditional Heteroscedasticity. \emph{arXiv:1609.00711}
}
\author{
Philipp Otto \email{potto@europa-uni.de}
}
\section{Control Arguments}{
\itemize{
  \item \code{seed} - {positive integer to initialize the random number generator (RNG), default value is a random integer in \eqn{[1, 10^6]}{[1, 10^6]}}
  \item \code{silent} - {if \code{FALSE}, current random seed is reported}
  \item \code{triangular} - {if \code{TRUE}, \eqn{\mathbf{W}}{W} is a triangular matrix and there are no checks to verify this assumption (default \code{FALSE})}
}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require("spdep")

# 1st example
##############

# parameters

rho <- 0.5
alpha <- 1
d <- 2

nblist <- cell2nb(d, d, type = "queen")
W <- nb2mat(nblist)

# simulation

Y <- sim.spARCH(rho = rho, alpha = alpha, W = W, type = "exp")

# visualization

image(1:d, 1:d, array(Y, dim = c(d,d)), xlab = expression(s[1]), ylab = expression(s[2]))

# 2nd example
##############

# two spatial weighting matrices W_1 and W_2
# h = alpha + rho_1 W_1 Y^2 + rho_2 W_2 Y^2

W_1 <- W
nblist <- cell2nb(d, d, type = "rook")
W_2 <- nb2mat(nblist)

rho_1 <- 0.3
rho_2 <- 0.7

W <- rho_1 * W_1 + rho_2 * W_2
rho <- 1

Y <- sim.spARCH(n = d^2, rho = rho, alpha = alpha, W = W, type = "exp")
image(1:d, 1:d, array(Y, dim = c(d,d)), xlab = expression(s[1]), ylab = expression(s[2]))
}
\keyword{spatial ARCH}
\keyword{spARCH}
