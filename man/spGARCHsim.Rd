\name{sim.spGARCH}
\alias{sim.spGARCH}
\title{Simulation of spatial ARCH models}
\description{
The function generates \code{n} random numbers of a spatial GARCH process for given parameters and weighting schemes.
}
\usage{
sim.spGARCH(n = dim(W1)[1], rho, lambda, alpha, W1, W2,
            b = 2, zeta = 0.5, theta = 0.5, type = "spGARCH",
            control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of observations. If \code{length(n)} \eqn{> 1}{> 1}, the length is taken to be the number required. Default \code{dim(W1)[1]}}
  \item{rho}{spatial dependence parameter rho}
  \item{lambda}{spatial dependence parameter lambda}
  \item{alpha}{unconditional variance level alpha}
  \item{W1}{\code{n} times \code{n} spatial weight matrix (ARCH component, parameter rho)}
  \item{W2}{\code{n} times \code{n} spatial weight matrix (GARCH component, parameter lambda)}
  \item{b}{parameter \code{b} for logarithmic spatial GARCH (only needed if \code{type = "log-spGARCH"}). Default 2.}
  \item{zeta}{parameter \code{zeta} for exponential spatial GARCH (only needed if \code{type = "e-spGARCH"}). Default 0.5.}
  \item{theta}{parameter \code{theta} for exponential spatial GARCH (only needed if \code{type = "e-spGARCH"}). Default 0.5.}
  \item{type}{type of simulated spGARCH process (see details)}
  \item{control}{list of control arguments (see below)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function simulates \code{n} observations \eqn{Y = (Y_1, ..., Y_n)'}{Y = (Y_1, ..., Y_n)'} of a spatial GARCH process, i.e.,
\deqn{
\boldsymbol{Y} = diag(\boldsymbol{h})^{1/2} \boldsymbol{\varepsilon} \, ,
}{
Y = diag(h)^(1/2) \epsilon ,
}
where \eqn{\boldsymbol{\varepsilon}}{\epsilon} is a spatial White Noise process. The definition of \eqn{\boldsymbol{h}}{h} depends on the chosen \code{type}. The following types are available.
\itemize{
\item \code{type = "spGARCH"} - simulates \eqn{\boldsymbol{\varepsilon}}{\epsilon} from a truncated normal distribution on the interval \eqn{[-a, a]}{[-a, a]}, such that \eqn{\boldsymbol{h} > 0}{h > 0} with
\deqn{
\boldsymbol{h} =  \alpha + \rho \mathbf{W}_1 \boldsymbol{Y}^{(2)} + \lambda \mathbf{W}_2 \boldsymbol{h} \; \mbox{and} \; a = 1 / ||\rho^2\mathbf{W}_1^2||_1^{1/4}.
}{
h = \alpha + \rho W1 Y^(2) + \lambda W2 h  and  a = 1 / (\rho^2||W1^2||_1)^(1/4)
}
Note that the normal distribution is not trunctated (\eqn{a = \infty}{a =  \infty}), if \eqn{\mathbf{W}_1}{W1} is a strictly triangular matrix, as it is ensured that \eqn{\boldsymbol{h} > \boldsymbol{0}}{h > 0}. Generally, it is sufficient that if there exists a permutation such that \eqn{\mathbf{W}_1}{W1} is strictly triangular. In this case, the process is called oriented spGARCH process.
\item \code{type = "e-spGARCH"} - simulates an exponential spARCH process (e-spGARCH), i.e.,
\deqn{
\ln\boldsymbol{h} = \alpha + \rho \mathbf{W}_1 g(\boldsymbol{\varepsilon}) + \lambda \mathbf{W}_2 log(\boldsymbol{h}) \, .
}{
ln h = \alpha + \rho W1 g(\epsilon) + \lambda W2 log(h).
}
For the e-spGARCH process, the errors follow a standard normal distribution. The function \eqn{g}{g} is given by
\deqn{
g(\boldsymbol{\varepsilon}) = \Theta {\varepsilon} + \zeta (|\varepsilon| - E(|\varepsilon|))  \, .
}{
g(\epsilon) = \Theta \epsilon + \zeta (|\epsilon| - E(|\epsilon|)) .
}
\item \code{type = "log-spGARCH"} - simulates a logarithmic spARCH process (log-spGARCH), i.e.,
\deqn{
\ln\boldsymbol{h} = \alpha + \rho \mathbf{W}_1 g(\boldsymbol{\varepsilon}) + \lambda \mathbf{W}_2 log(\boldsymbol{h}) \, .
}{
ln h = \alpha + \rho W1 g_b(\epsilon) + \lambda W2 log(h).
}
For the log-spGARCH process, the errors follow a standard normal distribution. The function \eqn{g}{g} is given by
\deqn{
g(\boldsymbol{\varepsilon}) = (\ln|\varepsilon(\boldsymbol{s}_1)|^{b}, \ldots, \ln|\varepsilon(\boldsymbol{s}_n)|^{b})' \, .
}{
g(\epsilon) = (ln|h(\epsilon_1)|^b, ..., ln|h(\epsilon_n)|^b)' .
}
\item \code{type = "complex-spGARCH"} - allows for complex solutions of \eqn{\boldsymbol{h}^{1/2}}{h^(1/2)} with
\deqn{
\boldsymbol{h} =  \alpha + \rho \mathbf{W}_1 \boldsymbol{Y}^{(2)} + \lambda \mathbf{W}_2 \boldsymbol{h}  \, .
}{
h =  \alpha + \rho W1 Y^(2)  + \lambda W2 h .
}
The errors follow a standard normal distribution.
}
}
\value{
The functions returns a vector \eqn{\boldsymbol{y}}{y}.
}
\references{
Philipp Otto, Wolfgang Schmid (2019). Spatial GARCH Models - A Unified Approach. \emph{arXiv:1908.08320}
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

# 1st example (spatial GARCH)
##############

# parameters

rho    <- 0.5
lambda <- 0.3
alpha  <- 1
d      <- 5

nblist <- cell2nb(d, d, type = "rook") # lattice process with Rook's contiguity matrix
W_1    <- nb2mat(nblist)
W_2    <- W_1

# simulation

Y      <- sim.spGARCH(rho = rho, lambda = lambda, alpha = alpha,
                      W1 = W_1, W2 = W_2, type = "spGARCH")

# visualization

image(1:d, 1:d, array(Y, dim = c(d,d)), xlab = expression(s[1]), ylab = expression(s[2]))

# 2nd example (exponential spatial GARCH)
##############

# parameters

rho    <- 0.5
lambda <- 0.3
alpha  <- 1
zeta   <- 0.5
theta  <- 0.5
d      <- 5

nblist <- cell2nb(d, d, type = "rook") # lattice process with Rook's contiguity matrix
W_1    <- nb2mat(nblist)
W_2    <- W_1

# simulation

Y      <- sim.spGARCH(rho = rho, lambda = lambda, alpha = alpha,
                      W1 = W_1, W2 = W_2, zeta = zeta, theta = 0.5, type = "e-spGARCH")

# visualization

image(1:d, 1:d, array(Y, dim = c(d,d)), xlab = expression(s[1]), ylab = expression(s[2]))

}
\keyword{spatial GARCH}
\keyword{spGARCH}
