\name{qml.spARCH}
\alias{qml.spARCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum-likelihood estimation of a spatial ARCH model
}
\description{
The function fits a spatial ARCH model using the maximum-likelihood approach. In addition, external regressor may be included in the mean equation.
}
\usage{
qml.spARCH(formula, W, type = "spARCH", data = NULL,
            b = 2, start = NULL, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class "\code{\link[stats]{formula}}"): a symbolic description of the model to be fitted. The details of model specification are given under "Details.""  \code{y} must be a numeric vector of length \code{n} (dimension of the weighting matrix)
}
  \item{W}{
\code{n} times \code{n} spatial weight matrix
}
  \item{type}{
type of spatial ARCH model to be fitted (see Details)
}
  \item{data}{
an optional data frame, list or environment containing the variables in the model. If not found in data, the variables are taken from the working space.
}
  \item{b}{
parameter \code{b} for the E-spARCH model (\code{type = "exp"}), must be a (positive) natural number
}
  \item{start}{
vector of starting values for the numerical optimization of the log-likelihood (optional)
}
  \item{control}{
list of control variables for iterative maximization of the log-likelihood
}
}
\details{
For \code{type = "spARCH"}, the functions fits a simple spatial ARCH model with one spatial lag, i.e.,
\deqn{
\boldsymbol{Y} = \mathbf{X} \boldsymbol{\beta} + \boldsymbol{h}^{1/2} \boldsymbol{\varepsilon}
}{
Y = X \beta + h^(1/2) \epsilon
}
with
\deqn{
\boldsymbol{h} = \alpha \boldsymbol{1}_n + \rho \mathbf{W} \boldsymbol{Y}^{(2)} \, .
}{
h = \alpha + \rho * W * Y^2 .
}
The distribution of the error term is assumed to be spARCH.

If \code{type = "exp"}, a spatial E-ARCH process is estimated, i.e.,
\deqn{
\ln(\boldsymbol{h}) = \alpha \boldsymbol{1}_n + \rho \mathbf{W} g_b(\boldsymbol{\varepsilon}) \, .
}{
ln h = \alpha + \rho * W * g_b(\epsilon) .
}
The function \eqn{g_b}{g_b} is defined as
\deqn{
g_b(\boldsymbol{\varepsilon}) = (\ln|\varepsilon(\boldsymbol{s}_1)|^{b}, \ldots, \ln|\varepsilon(\boldsymbol{s}_n)|^{b})'
}{
g_b(\epsilon) = (ln|h(\epsilon_1)|^b, ..., ln|h(\epsilon_n)|^b)'
}
and the error term is also assumed to be spARCH.

The modelling equation gan be specified as for \code{\link{lm}}, i.e., as \code{\link[stats]{formula}} object. A typical model has the form \code{response ~ terms} where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for \code{response}. A terms specification of the form \code{first + second} indicates all the terms in first together with all the terms in second with duplicates removed. A specification of the form \code{first:second} indicates the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification \code{first*second} indicates the cross of first and second. This is the same as \code{first + second + first:second}. However, there is no \code{offset} permitted for the \code{\link{qml.spARCH}}.

For an intercept-only model, the \code{\link[stats]{formula}} can be specified as \code{response ~ 1}. In addition, it is possible to fit an intercept-free model with \code{response ~ 0} or \code{response ~ 0 + terms}.

To summarize the results of the model fit, use the generic function \code{summary}. For analysis of the residuals, the generic \code{plot} provides several descriptive plots. For numerical maximization of the log-likelihood, the function uses the algorithm of \code{\link[Rsolnp]{solnp}} from the package \code{\link[Rsolnp]{Rsolnp}}.
}
\value{
A \linkS4class{spARCH} object with the following elements is returned:
\item{coefficients}{Parameter estimates \eqn{\alpha}{\alpha} and \eqn{\rho}{\rho}.}
\item{residuals}{Vector of residuals.}
\item{fitted.values}{Fitted values.}
\item{stderr}{Standard errors of the estimates (Cramer-Rao estimates).}
\item{hessian}{Hessian matrix of the negative Log-Likelihood at the estimated minimum.}
\item{LL}{Value of the Log-Likelihood at the estimated maximum.}
\item{h}{Fitted vector \eqn{\boldsymbol{h}}{h}.}
\item{y}{Vector of observations (input values).}
\item{h}{Chosen type (input).}
\item{W}{Spatial weight matrix (input).}
\item{regressors}{Are regressors included? \code{TRUE}/\code{FALSE}}
\item{AR}{Is an autoregressive term in the mean equation? \code{TRUE}/\code{FALSE}}
\item{X}{Matrix of regressors if \code{regressor = TRUE}}
}
\references{
Philipp Otto, Wolfgang Schmid, Robert Garthoff (2018). Generalised Spatial and Spatiotemporal Autoregressive Conditional Heteroscedasticity. \href{https://doi.org/10.1016/j.spasta.2018.07.005}{Spatial Statistics} 26, pp. 125-145. \emph{arXiv:1609.00711}
}
\author{
Philipp Otto \email{potto@europa-uni.de}
}
\section{Control Arguments}{
see also: \code{\link{solnp}}
\itemize{
  \item \code{rho} - {This is used as a penalty weighting scalar for infeasibility in the augmented objective function. The higher its value the more the weighting to bring the solution into the feasible region (default 1). However, very high values might lead to numerical ill conditioning or significantly slow down convergence.}
  \item \code{outer.iter} - {Maximum number of major (outer) iterations (default 400).}
  \item \code{inner.iter} - {Maximum number of minor (inner) iterations (default 800).}
  \item \code{delta} - {Relative step size in forward difference evaluation (default 1.0e-7).}
  \item \code{tol} - {Relative tolerance on feasibility and optimality (default 1e-8).}
  \item \code{trace} - {The value of the objective function and the parameters is printed at every major iteration (default 1).}
}
}
\seealso{
\code{\link{solnp}} (package \code{\link[Rsolnp]{Rsolnp}})
}
\examples{
require("spdep")

# directional spatial ARCH process (W is triangular, 1:1 origin)

rho <- 0.5
alpha <- 1
d <- 5
n <- d^2
nblist <- cell2nb(d, d, type = "queen")
W <- nb2mat(nblist)
W[lower.tri(W)] <- 0
y <- sim.spARCH(n = n, rho = rho, alpha = alpha, W = W, type = "spARCH")

out <- qml.spARCH(y ~ 0, W = W, type = "spARCH")

summary(out)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Spatial ARCH Model}% use one of  RShowDoc("KEYWORDS")
\keyword{QML spARCH}% __ONLY ONE__ keyword per line
