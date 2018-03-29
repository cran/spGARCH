#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
const double pi = 3.14159265358979323846;


//-----------------------------//
//                             //
//    spatial ARCH models      //
//                             //
//-----------------------------//



//-----------------------------//
//                             //
//    Models w/o regressors    //
//                             //
//-----------------------------//


// [[Rcpp::export(".LL_spARCH")]]
double LL_spARCH(Eigen::VectorXd pars, List param) {

  // ---------------- //
  //   Declarations   //
  // ---------------- //


  Eigen::VectorXd                         y = param[0];
  Eigen::MappedSparseMatrix<double> Wsparse = param[1];
  bool                                  sym = param[2];
  int                                     n = y.size();
  double                              alpha = pars[0];
  double                                rho = pars[1];

  // ---------------- //
  //   Computation    //
  // ---------------- //


  VectorXd                h = alpha * VectorXd::Ones(n) + rho * Wsparse * y.asDiagonal() * y;
  ArrayXd               eps = y.array() / (h.array()).sqrt();

  ArrayXd               aux = h.array() / pow(y.array(), 2);
  MatrixXd             aux2 = (aux.matrix()).asDiagonal();
  MatrixXd             rhoW = rho * Wsparse;

  VectorXd eigen_v(n);
  if(sym)
  {
    SelfAdjointEigenSolver<MatrixXd>       es(aux2 - rhoW.transpose());
                                 eigen_v = es.eigenvalues();
  }
  else
  {
    EigenSolver<MatrixXd>       es(aux2 - rhoW.transpose());
                      eigen_v = es.eigenvalues().real();
  }

  double log_abs_det_partial = ( (pow(y.array(), 2) / (pow(h.array(), 3)).sqrt()).log() + ((eigen_v.array()).abs()).log()).sum();
  double sum_f_residuals     = - 0.5 * log(2 * pi) - 0.5 * (pow(eps.array(), 2)).sum();

  return  (-1) * (sum_f_residuals + log_abs_det_partial);
}

// [[Rcpp::export(".LL_EspARCH")]]
double LL_EspARCH(Eigen::VectorXd pars, List param) {

  // ---------------- //
  //   Declarations   //
  // ---------------- //

  VectorXd                                y = param[0];
  Eigen::MappedSparseMatrix<double> Wsparse = param[1];
  int                                     b = param[2];

  int                                     n = y.size();
  double                              alpha = pars[0];
  double                                rho = pars[1];
  MatrixXd                             rhoW = rho * Wsparse;

  // Identity matrix
  MatrixXd I(n, n);
  I.setIdentity(n, n);

  // ---------------- //
  //   Computation    //
  // ---------------- //

  VectorXd lnabsY = ((y.array()).abs()).log();

  MatrixXd S     = I + 0.5 * b * rhoW;
           S     = S.inverse();
  VectorXd foo   = alpha * VectorXd::Ones(n) + b * rhoW * lnabsY;
  VectorXd lnh   = S * foo;
  ArrayXd  h     = lnh.array();
           h     = h.exp();
  ArrayXd  eps   = y.array() / h.sqrt();
  ArrayXd  aux   = 2 * h / b;
  MatrixXd aux_2 = (aux.matrix()).asDiagonal();
  MatrixXd FOO   = (rhoW.array() * S.array()).matrix();

  EigenSolver<MatrixXd> es(aux_2 - FOO.transpose());
  VectorXd eigen_values = es.eigenvalues().real();
  ArrayXd       eigen_v = (eigen_values.array()).abs();

  double log_abs_det_partial = (log(b / (2 * sqrt(pow(h.array(), 3)))) + log(eigen_v)).sum();
  double sum_f_residuals     = - 0.5 * log(2 * pi) - 0.5 * (pow(eps.array(), 2)).sum();

  return  (-1) * (sum_f_residuals + log_abs_det_partial);
}

//-----------------------------//
//                             //
//     Models w/ regressors    //
//                             //
//-----------------------------//


// [[Rcpp::export(".LL_spARCHX")]]
double LL_spARCHX(Eigen::VectorXd pars, List param) {

  // ---------------- //
  //   Declarations   //
  // ---------------- //

  Eigen::VectorXd                         y = param[0];
  Eigen::MappedSparseMatrix<double> Wsparse = param[1];
  Eigen::MatrixXd                         X = param[2];
  bool                                  sym = param[3];
  int                                     n = y.size();
  int                                     p = X.cols();
  double                              alpha = pars[0];
  double                                rho = pars[1];

  // ---------------- //
  //   Computation    //
  // ---------------- //

  VectorXd beta(p);
  for(int i = 0; i < p; i++)
  {
    beta[i] = pars[2 + i];
  }
  VectorXd xi    = y - X * beta;
  VectorXd h     = alpha * VectorXd::Ones(n) + rho * Wsparse * xi.asDiagonal() * xi;
  VectorXd eps   = (xi.array() / sqrt(h.array())).matrix();
  ArrayXd  aux   = h.array() / pow(xi.array(), 2);
  MatrixXd aux_2 = (aux.matrix()).asDiagonal();
  MatrixXd rhoW  = rho * Wsparse;

  ArrayXd  eigen_v(n);
  if(sym)
  {
    SelfAdjointEigenSolver<MatrixXd> es(aux_2 - rhoW.transpose());
    VectorXd          eigen_values = es.eigenvalues();
                      eigen_v      = (eigen_values.array()).abs();
  }
  else
  {
    EigenSolver<MatrixXd> es(aux_2 - rhoW.transpose());
    VectorXd eigen_values = es.eigenvalues().real();
             eigen_v      = (eigen_values.array()).abs();
  }

  double log_abs_det_partial = (log(pow(xi.array(), 2) / sqrt(pow(h.array(), 3))) + log(eigen_v)).sum();
  double sum_f_residuals     = - 0.5 * log(2 * pi) - 0.5 * (pow(eps.array(), 2)).sum();

  return  (-1) * (sum_f_residuals + log_abs_det_partial);
}

// [[Rcpp::export(.LL_EspARCHX)]]
double LL_EspARCHX(Eigen::VectorXd pars, List param) {

  // ---------------- //
  //   Declarations   //
  // ---------------- //

  VectorXd                                y = param[0];
  Eigen::MappedSparseMatrix<double> Wsparse = param[1];
  int                                     b = param[2];
  Eigen::MatrixXd                         X = param[3];

  int                                     n = y.size();
  int                                     p = X.cols();
  double                              alpha = pars[0];
  double                                rho = pars[1];

  VectorXd lnabsxi(n);
  MatrixXd rhoW = rho * Wsparse;

  // Identity matrix
  MatrixXd I(n, n);
  I.setIdentity(n, n);

  VectorXd beta(p);
  for(int i = 0; i < p; i++)
  {
    beta[i] = pars[2 + i];
  }


  // ---------------- //
  //   Computation    //
  // ---------------- //

  VectorXd xi  = y - X * beta;
  lnabsxi      = ((xi.array()).abs()).log();

  MatrixXd S   = I + 0.5 * b * rhoW;
           S   = S.inverse();
  VectorXd foo = alpha * VectorXd::Ones(n) + b * rhoW * lnabsxi;
  VectorXd lnh = S * foo;
  ArrayXd    h = lnh.array();
             h = h.exp();
  ArrayXd  eps   = xi.array() / h.sqrt();
  ArrayXd  aux   = 2 * h / b;
  MatrixXd aux_2 = (aux.matrix()).asDiagonal();
  MatrixXd FOO   = (rhoW.array() * S.array()).matrix();

  EigenSolver<MatrixXd>   es(aux_2 - FOO.transpose());
  VectorXd eigen_values = es.eigenvalues().real();
  ArrayXd  eigen_v      = (eigen_values.array()).abs();

  double log_abs_det_partial = (log(b / (2 * sqrt(pow(h.array(), 3)))) + log(eigen_v)).sum();
  double sum_f_residuals     = - 0.5 * log(2 * pi) - 0.5 * (pow(eps.array(), 2)).sum();

  return  (-1) * (sum_f_residuals + log_abs_det_partial);
}






//-----------------------------//
//                             //
//      SARspARCH models       //
//                             //
//-----------------------------//

//-----------------------------//
//                             //
//    Models w/o regressors    //
//                             //
//-----------------------------//


// [[Rcpp::export(".LL_SARspARCH")]]
double LL_SARspARCH(Eigen::VectorXd pars, List param) {

  // ---------------- //
  //   Declarations   //
  // ---------------- //


  Eigen::VectorXd                         y = param[0];
  Eigen::MappedSparseMatrix<double> Bsparse = param[1];
  Eigen::MappedSparseMatrix<double> Wsparse = param[2];
  Eigen::VectorXd                   B_eigen = param[3];
  bool                                  sym = param[4];
  int                                     n = y.size();
  double                              alpha = pars[0];
  double                                rho = pars[1];
  double                             lambda = pars[2];

  // ---------------- //
  //   Computation    //
  // ---------------- //

  VectorXd xi    = y - lambda * Bsparse * y;
  VectorXd h     = alpha * VectorXd::Ones(n) + rho * Wsparse * xi.asDiagonal() * xi;
  VectorXd eps   = (xi.array() / sqrt(h.array())).matrix();
  ArrayXd  aux   = h.array() / pow(xi.array(), 2);
  MatrixXd aux_2 = (aux.matrix()).asDiagonal();
  MatrixXd rhoW  = rho * Wsparse;

  ArrayXd  eigen_v(n);
  if(sym)
  {
    SelfAdjointEigenSolver<MatrixXd> es(aux_2 - rhoW.transpose());
    VectorXd          eigen_values = es.eigenvalues();
                      eigen_v      = (eigen_values.array()).abs();
  }
  else
  {
    EigenSolver<MatrixXd> es(aux_2 - rhoW.transpose());
    VectorXd eigen_values = es.eigenvalues().real();
             eigen_v      = (eigen_values.array()).abs();
  }

  double log_abs_det_partial    = (log(pow(xi.array(), 2) / sqrt(pow(h.array(), 3))) + log(eigen_v)).sum();
  double log_abs_det_partial_ar = (((1 - lambda * B_eigen.array()).abs()).log()).sum();
  double sum_f_residuals        = - 0.5 * log(2 * pi) - 0.5 * (pow(eps.array(), 2)).sum();

  return  (-1) * (sum_f_residuals + log_abs_det_partial + log_abs_det_partial_ar);
}

// [[Rcpp::export(".LL_SAREspARCH")]]
double LL_SAREspARCH(Eigen::VectorXd pars, List param) {

  // ---------------- //
  //   Declarations   //
  // ---------------- //

  VectorXd                                y = param[0];
  Eigen::MappedSparseMatrix<double> Bsparse = param[1];
  Eigen::MappedSparseMatrix<double> Wsparse = param[2];
  int                                     b = param[3];
  Eigen::VectorXd                   B_eigen = param[4];

  int                                     n = y.size();
  double                              alpha = pars[0];
  double                                rho = pars[1];
  double                             lambda = pars[2];

  VectorXd lnabsxi(n);
  MatrixXd                             rhoW = rho * Wsparse;

  // Identity matrix
  MatrixXd I(n, n);
  I.setIdentity(n, n);

  // ---------------- //
  //   Computation    //
  // ---------------- //

  VectorXd xi = y - lambda * Bsparse * y;
      lnabsxi = ((xi.array()).abs()).log();

  MatrixXd S   = I + 0.5 * b * rhoW;
           S   = S.inverse();
  VectorXd foo = alpha * VectorXd::Ones(n) + b * rhoW * lnabsxi;
  VectorXd lnh = S * foo;
  ArrayXd    h = lnh.array();
             h = h.exp();
  ArrayXd  eps   = xi.array() / h.sqrt();
  ArrayXd  aux   = 2 * h / b;
  MatrixXd aux_2 = (aux.matrix()).asDiagonal();
  MatrixXd FOO   = (rhoW.array() * S.array()).matrix();

  EigenSolver<MatrixXd>   es(aux_2 - FOO.transpose());
  VectorXd eigen_values = es.eigenvalues().real();
  ArrayXd  eigen_v      = (eigen_values.array()).abs();

  double log_abs_det_partial    = (log(b / (2 * sqrt(pow(h.array(), 3)))) + log(eigen_v)).sum();
  double log_abs_det_partial_ar = (((1 - lambda * B_eigen.array()).abs()).log()).sum();
  double sum_f_residuals        = - 0.5 * log(2 * pi) - 0.5 * (pow(eps.array(), 2)).sum();

  return  (-1) * (sum_f_residuals + log_abs_det_partial + log_abs_det_partial_ar);
}



//-----------------------------//
//                             //
//     Models w/ regressors    //
//                             //
//-----------------------------//



// [[Rcpp::export(".LL_SARspARCHX")]]
double LL_SARspARCHX(Eigen::VectorXd pars, List param) {

  // ---------------- //
  //   Declarations   //
  // ---------------- //

  Eigen::VectorXd                         y = param[0];
  Eigen::MappedSparseMatrix<double> Bsparse = param[1];
  Eigen::MappedSparseMatrix<double> Wsparse = param[2];
  Eigen::MatrixXd                         X = param[3];
  Eigen::VectorXd                   B_eigen = param[4];
  bool                                  sym = param[5];

  double                              alpha = pars[0];
  double                                rho = pars[1];
  double                             lambda = pars[2];

  int n = y.size();
  int p = X.cols();

  VectorXd beta(p);
  for(int i = 0; i < p; i++)
  {
    beta[i] = pars[3 + i];
  }

  // ---------------- //
  //   Computation    //
  // ---------------- //


  VectorXd xi    = y - lambda * Bsparse * y - X * beta;
  VectorXd h     = alpha * VectorXd::Ones(n) + rho * Wsparse * xi.asDiagonal() * xi;
  VectorXd eps   = (xi.array() / sqrt(h.array())).matrix();
  ArrayXd  aux   = h.array() / pow(xi.array(), 2);
  MatrixXd aux_2 = (aux.matrix()).asDiagonal();
  MatrixXd rhoW  = rho * Wsparse;

  ArrayXd  eigen_v(n);
  if(sym)
  {
    SelfAdjointEigenSolver<MatrixXd> es(aux_2 - rhoW.transpose());
    VectorXd          eigen_values = es.eigenvalues();
                      eigen_v      = eigen_values.array().abs();
  }
  else
  {
    EigenSolver<MatrixXd> es(aux_2 - rhoW.transpose());
    VectorXd eigen_values = es.eigenvalues().real();
             eigen_v      = eigen_values.array().abs();
  }

  double log_abs_det_partial    = (log(pow(xi.array(), 2) / sqrt(pow(h.array(), 3))) + log(eigen_v)).sum();
  double log_abs_det_partial_ar = (((1 - lambda * B_eigen.array()).abs()).log()).sum();
  double sum_f_residuals        = - 0.5 * log(2 * pi) - 0.5 * (pow(eps.array(), 2)).sum();

  return  (-1) * (sum_f_residuals + log_abs_det_partial + log_abs_det_partial_ar);
}

// [[Rcpp::export(".LL_SAREspARCHX")]]
double LL_SAREspARCHX(Eigen::VectorXd pars, List param) {

  // ---------------- //
  //   Declarations   //
  // ---------------- //

  VectorXd                                y = param[0];
  Eigen::MappedSparseMatrix<double> Bsparse = param[1];
  Eigen::MappedSparseMatrix<double> Wsparse = param[2];
  int                                     b = param[3];
  Eigen::MatrixXd                         X = param[4];
  Eigen::VectorXd                   B_eigen = param[5];

  double                              alpha = pars[0];
  double                                rho = pars[1];
  double                             lambda = pars[2];

  int n = y.size();
  int p = X.cols();

  VectorXd lnabsxi(n);
  MatrixXd rhoW = rho * Wsparse;

  // Identity matrix
  MatrixXd I(n, n);
  I.setIdentity(n, n);

  VectorXd beta(p);
  for(int i = 0; i < p; i++)
  {
    beta[i] = pars[3 + i];
  }


  // ---------------- //
  //   Computation    //
  // ---------------- //

  VectorXd       xi = y - lambda * Bsparse * y - X * beta;
            lnabsxi = ((xi.array()).abs()).log();

  MatrixXd   S = (I + 0.5 * b * rhoW).inverse();
  VectorXd foo = alpha * VectorXd::Ones(n) + b * rhoW * lnabsxi;
  VectorXd lnh = S * foo;
  ArrayXd    h = (lnh.array()).exp();

  ArrayXd  eps   = xi.array() / h.sqrt();
  ArrayXd  aux   = 2 * h / b;
  MatrixXd aux_2 = (aux.matrix()).asDiagonal();
  MatrixXd FOO   = (rhoW.array() * S.array()).matrix();

  EigenSolver<MatrixXd>   es(aux_2 - FOO.transpose());
  VectorXd eigen_values = es.eigenvalues().real();
  ArrayXd  eigen_v      = (eigen_values.array()).abs();

  double log_abs_det_partial    = (log(b / (2 * sqrt(pow(h.array(), 3)))) + log(eigen_v)).sum();
  double log_abs_det_partial_ar = (((1 - lambda * B_eigen.array()).abs()).log()).sum();
  double sum_f_residuals        = - 0.5 * log(2 * pi) - 0.5 * (pow(eps.array(), 2)).sum();

  return  (-1) * (sum_f_residuals + log_abs_det_partial + log_abs_det_partial_ar);
}
