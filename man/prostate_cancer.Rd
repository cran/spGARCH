\name{prostate_cancer}
\alias{prostate_cancer}
\docType{data}
\title{
Logarithmic incidence rates of prostate cancer and covariates
}
\description{
The dataset contains logarithmic incidence rates from the National Cancer Institute and Centers for Disease, Control and Prevention (State Cancer Profiles) and factor loadings of a set of covariates. The incidence rates are 5-year averages from 2008 to 2012 in several southeastern states (Arkansas, Louisiana, Mississippi, Tennessee, North and South Carolina, Georgia, Alabama, and Florida). Missing values were imputed by spatial averaging.
}
\usage{data("prostate_cancer")}
\format{
  A list with three entries:
  \describe{
    \item{\code{data}}{a data frame; see below for details}
    \item{\code{B}}{a numeric matrix; weighting matrix B to run the example}
    \item{\code{W}}{a numeric matrix; weighting matrix W to run the example}
  }
  The data frame contains 755 observations of the following 12 variables.
  \describe{
    \item{\code{log_incidence_rates}}{a numeric vector; logarithmic incidence rates of prostate cancer}
    \item{\code{F_1}}{a numeric vector; scores of factor 1 (environment: fine atmospheric particles and aerosols)}
    \item{\code{F_2}}{a numeric vector; scores of factor 2 (environment: particulate matter)}
    \item{\code{F_3}}{a numeric vector; scores of factor 3 (weather: solar radiation and temperature)}
    \item{\code{F_4}}{a numeric vector; scores of factor 4 (weather: temperature differences)}
    \item{\code{F_5}}{a numeric vector; scores of factor 5 (behavior: smoking)}
    \item{\code{F_6}}{a numeric vector; scores of factor 6 (behavior: drinking)}
    \item{\code{F_7}}{a numeric vector; scores of factor 7 (behavior: preventive health care)}
    \item{\code{F_8}}{a numeric vector; scores of factor 8 (behavior: physical activity)}
    \item{\code{F_9}}{a numeric vector; scores of factor 9 (health: overweight)}
    \item{\code{F_10}}{a numeric vector; scores of factor 10 (health: cholesterol and blood pressure)}
    \item{\code{PSA_test}}{a numeric vector;  percentage of positive results for a prostate-specific antigen (PSA) test}
  }
}
\source{
\url{https://statecancerprofiles.cancer.gov/index.html}

National Cancer Institute, Centers for Disease, Control and Prevention
}
\references{
https://statecancerprofiles.cancer.gov/map/map.withimage.php?99&001&001&00&0&02&0&1&10
Otto, P. (2019). spGARCH: An R-Package for Spatial and Spatiotemporal ARCH and GARCH models To appear: The R Journal URL: https://arxiv.org/abs/1812.01871
}
\examples{
data(prostate_cancer)

\dontrun{

data <- prostate_cancer$data

# Estimation (long running example)
formula <- "log_incidence_rates ~ F_2 + F_10"
out     <- qml.SARspARCH(formula, prostate_cancer$B = B, W = prostate_cancer$W,
                         type = "gaussian", data = data)

# Summary
summary(out)

}

}
\keyword{datasets}
