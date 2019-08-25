.onAttach <- function(...) {
  packageStartupMessage(
"\n Generalized spatial autoregressive conditional heteroscedasticity (spGARCH) \n
 Theoretical details of spatial ARCH models can be found in Otto, Schmid, Garthoff (2018): \n
\t Otto, P., Schmid, W. & Garthoff, R. (2018)
\t Generalised spatial and spatiotemporal autore-
\t gressive conditional heteroscedasticity,
\t Spatial Statistics 26, pp. 125-145,
\t arXiv:1609.00711 \n
Regarding spatial GARCH models, we refer to Otto, Schmid (2019): \n
\t Otto, P. & Schmid, W. (2019)
\t Spatial GARCH models - A unified approach
\t arXiv:1908.08320 \n
 Computational implementation: \n
\t Otto, P. (2019)
\t spGARCH: An R-package for spatial and spatiotemporal ARCH
\t and GARCH models
\t To appear in: The R Journal,
\t arXiv:1812.01871v1 \n")
}
# To suppress this message use:
#   suppressPackageStartupMessages(library(installr))



