.onAttach <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    blue <- crayon::blue
    green <- crayon::green
    bold <- crayon::bold
  } else {
    blue <- green <- bold <- identity
  }
  
  packageStartupMessage(
    bold("\n--- spGARCH: Generalized Spatial ARCH and GARCH Models ---\n"),
    "\n", green("Theoretical details of spatial ARCH models:"), "\n",
    "  ", blue("Otto, P., Schmid, W., Garthoff, R. (2018)."), 
    " Generalised spatial and spatiotemporal\n",
    "  autoregressive conditional heteroscedasticity.\n", 
    "  Spatial Statistics, 26, 125-145.\n",
    "  DOI: 10.1016/j.spasta.2018.07.005 | arXiv:1609.00711\n",
    
    "\n", green("Review paper on spatial and spatiotemporal volatility models:"), "\n",
    "  ", blue("Otto, P., Dogan, O., Taspinar, S., Schmid, W., Bera, A. K. (2019)."), 
    "  Spatial and spatiotemporal volatility models: A review.\n",
    "  Journal of Economic Surveys\n",
    "  DOI: 10.1111/joes.12643 | arXiv:2308.13061\n",
    
    "\n", green("For implementation details:"), "\n",
    "  ", blue("Otto, P. (2019)."), 
    "  spGARCH: An R Package for Spatial and Spatiotemporal ARCH and GARCH Models.\n",
    "  The R Journal, 11(2).\n",
    "  URL: https://journal.r-project.org/articles/RJ-2019-053/ | arXiv:1812.01871 \n"
  )
}