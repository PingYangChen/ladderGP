
#include "header.h"

// RCPP FUNCTIONS

//[[Rcpp::export]]
double gpObjCpp(arma::vec param, arma::vec y, arma::mat x, double nugget)
{
  arma::uword n = x.n_rows;
  double negloglik, mu, sigma2;
  arma::mat psi(n, n, fill::eye);
  arma::mat invPsi(n, n, fill::eye);
  //double nugget = 0.;
  gpLogLik(negloglik, psi, invPsi, mu, sigma2, nugget, y, x, param);
  return negloglik;
}

//[[Rcpp::export]]
Rcpp::List gpModel(arma::vec param, arma::vec y, arma::mat x, double nugget)
{
  arma::uword n = x.n_rows;
  double negloglik, mu, sigma2;
  arma::mat psi(n, n, fill::eye);
  arma::mat invPsi(n, n, fill::eye);
  //double nugget = 0.;
  gpLogLik(negloglik, psi, invPsi, mu, sigma2, nugget, y, x, param);
  return List::create(Named("alpha") = wrap(param),
                      Named("psi") = wrap(psi),
                      Named("invPsi") = wrap(invPsi),
                      Named("mu") = wrap(mu),
                      Named("sigma2") = wrap(sigma2),
                      Named("negloglik") = wrap(negloglik),
                      Named("nugget") = wrap(nugget));
}

//[[Rcpp::export]]
Rcpp::List gpPred(arma::mat x0, arma::vec y, arma::mat x, 
                  arma::vec param, arma::mat invPsi, double mu, double sigma2, double ei_alpha, double min_y)
{
  arma::uword n0 = x0.n_rows;
  arma::vec y0(n0, fill::zeros);
  arma::vec mse(n0, fill::zeros);
  arma::vec ei(n0, fill::zeros);  
  arma::vec ei_1(n0, fill::zeros);
  arma::vec ei_2(n0, fill::zeros);
  gpNewData(y0, mse, ei, ei_1, ei_2, ei_alpha, min_y, x0, y, x, mu, sigma2, invPsi, param);
  //
  return List::create(Named("pred") = wrap(y0),
                      Named("mse") = wrap(mse),
                      Named("ei") = wrap(ei),
                      Named("improvement") = wrap(ei_1),
                      Named("uncertainty") = wrap(ei_2)
                     );
}


// Stage as a ordinal variable
//[[Rcpp::export]]
double lgpOBObjCpp(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget, bool logParam)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu, sigma2;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  double tau;
  arma::mat theta;
  //double nugget = 0.;
  if (logParam) {
    lgpOBParam2vec(tau, theta, arma::exp(param), xzDim, zMax);
  } else{
    lgpOBParam2vec(tau, theta, param, xzDim, zMax);
  }
  lgpOBLogLik(negloglik, psi, invPsi, mu, sigma2, nugget, y, x, z, xzDim, 
              tau, theta);
  return negloglik;
}

//[[Rcpp::export]]
Rcpp::List lgpOBModel(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget, bool logParam)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu, sigma2;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  double tau;
  arma::mat theta;
  //double nugget = 0.;
  if (logParam) {
    lgpOBParam2vec(tau, theta, arma::exp(param), xzDim, zMax);
  } else{
    lgpOBParam2vec(tau, theta, param, xzDim, zMax);
  }
  lgpOBLogLik(negloglik, psi, invPsi, mu, sigma2, nugget, y, x, z, xzDim, 
              tau, theta);
  //
  return List::create(Named("mu") = wrap(mu),
                      Named("sigma2") = wrap(sigma2),
                      Named("tau") = wrap(tau),
                      Named("theta") = wrap(theta),
                      Named("psi") = wrap(psi),
                      Named("invPsi") = wrap(invPsi),
                      Named("negloglik") = wrap(negloglik),
                      Named("nugget") = wrap(nugget),
                      Named("vecParams") = wrap(param)
  );
}

//[[Rcpp::export]]
Rcpp::List lgpOBPred(arma::mat x0, arma::uvec z0, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim,
                     arma::rowvec param, arma::mat invPsi, double mu, double sigma2, double ei_alpha, double min_y, bool logParam)
{
  arma::uword n0 = x0.n_rows;
  arma::uword zMax = z.max();
  //
  double tau;
  arma::mat theta;
  if (logParam) {
    lgpOBParam2vec(tau, theta, arma::exp(param), xzDim, zMax);
  } else{
    lgpOBParam2vec(tau, theta, param, xzDim, zMax);
  }
  //
  arma::vec y0(n0, fill::zeros);
  arma::vec mse(n0, fill::zeros);
  arma::vec ei(n0, fill::zeros);
  arma::vec ei_1(n0, fill::zeros);
  arma::vec ei_2(n0, fill::zeros);
  lgpOBNewData(y0, mse, ei, ei_1, ei_2, ei_alpha, min_y, x0, z0, y, x, z, xzDim, 
               mu, sigma2, invPsi, tau, theta);
  //
  return List::create(Named("pred") = wrap(y0),
                      Named("mse") = wrap(mse),
                      Named("ei") = wrap(ei),
                      Named("improvement") = wrap(ei_1),
                      Named("uncertainty") = wrap(ei_2)
  );
}


// Stage as a nominal variable

//[[Rcpp::export]]
double lgpNBObjCpp(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget, bool logParam)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu, sigma2;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  arma::mat tau, theta;
  //double nugget = 0.;
  if (logParam) {
    lgpNBParam2vec(tau, theta, arma::exp(param), xzDim, zMax);
  } else{
    lgpNBParam2vec(tau, theta, param, xzDim, zMax);
  }
  lgpNBLogLik(negloglik, psi, invPsi, mu, sigma2, nugget, y, x, z, xzDim, 
              tau, theta);
  return negloglik;
}

//[[Rcpp::export]]
Rcpp::List lgpNBModel(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget, bool logParam)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu, sigma2;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  arma::mat tau, theta;
  //double nugget = 0.;
  if (logParam) {
    lgpNBParam2vec(tau, theta, arma::exp(param), xzDim, zMax);
  } else{
    lgpNBParam2vec(tau, theta, param, xzDim, zMax);
  }
  lgpNBLogLik(negloglik, psi, invPsi, mu, sigma2, nugget, y, x, z, xzDim, 
              tau, theta);
  //
  return List::create(Named("mu") = wrap(mu),
                      Named("sigma2") = wrap(sigma2),
                      Named("tau") = wrap(tau),
                      Named("theta") = wrap(theta),
                      Named("psi") = wrap(psi),
                      Named("invPsi") = wrap(invPsi),
                      Named("negloglik") = wrap(negloglik),
                      Named("nugget") = wrap(nugget),
                      Named("vecParams") = wrap(param)
  );
}

//[[Rcpp::export]]
Rcpp::List lgpNBPred(arma::mat x0, arma::uvec z0, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim,
                     arma::rowvec param, arma::mat invPsi, double mu, double sigma2, double ei_alpha, double min_y, bool logParam)
{
  arma::uword n0 = x0.n_rows;
  arma::uword zMax = z.max();
  //
  arma::mat tau, theta;
  if (logParam) {
    lgpNBParam2vec(tau, theta, arma::exp(param), xzDim, zMax);
  } else{
    lgpNBParam2vec(tau, theta, param, xzDim, zMax);
  }
  //
  arma::vec y0(n0, fill::zeros);
  arma::vec mse(n0, fill::zeros);
  arma::vec ei(n0, fill::zeros);
  arma::vec ei_1(n0, fill::zeros);
  arma::vec ei_2(n0, fill::zeros);
  lgpNBNewData(y0, mse, ei, ei_1, ei_2, ei_alpha, min_y, x0, z0, y, x, z, xzDim, 
               mu, sigma2, invPsi, tau, theta);
  //
  return List::create(Named("pred") = wrap(y0),
                      Named("mse") = wrap(mse),
                      Named("ei") = wrap(ei),
                      Named("improvement") = wrap(ei_1),
                      Named("uncertainty") = wrap(ei_2)
  );
}


// Stage as interaction effect
//[[Rcpp::export]]
double aIntObjCpp(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget, bool logParam)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  arma::mat thetaZ;
  arma::vec sigmaF;
  arma::mat sigmaInt;
  //double nugget = 0.;
  if (logParam) {
    aIntParam2vec(thetaZ, sigmaF, sigmaInt, arma::exp(param), xzDim, zMax);
  } else{
    aIntParam2vec(thetaZ, sigmaF, sigmaInt, param, xzDim, zMax);
  }
  aIntLogLik(negloglik, psi, invPsi, mu, nugget, y, x, z, xzDim, 
             thetaZ, sigmaF, sigmaInt);
  return negloglik;
}

//[[Rcpp::export]]
Rcpp::List aIntModel(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget, bool logParam)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  arma::mat thetaZ;
  arma::vec sigmaF;
  arma::mat sigmaInt;
  //double nugget = 0.;
  if (logParam) {
    aIntParam2vec(thetaZ, sigmaF, sigmaInt, arma::exp(param), xzDim, zMax);
  } else{
    aIntParam2vec(thetaZ, sigmaF, sigmaInt, param, xzDim, zMax);
  }
  aIntLogLik(negloglik, psi, invPsi, mu, nugget, y, x, z, xzDim, 
             thetaZ, sigmaF, sigmaInt);
  //
  return List::create(Named("mu") = wrap(mu),
                      Named("thetaZ") = wrap(thetaZ),
                      Named("sigmaF") = wrap(sigmaF),
                      Named("sigmaInt") = wrap(sigmaInt),
                      Named("psi") = wrap(psi),
                      Named("invPsi") = wrap(invPsi),
                      Named("negloglik") = wrap(negloglik),
                      Named("nugget") = wrap(nugget),
                      Named("vecParams") = wrap(param)
  );
}

//[[Rcpp::export]]
Rcpp::List aIntPred(arma::mat x0, arma::uvec z0, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim,
                    arma::rowvec param, arma::mat invPsi, double mu, double ei_alpha, double min_y, bool logParam)
{
  arma::uword n0 = x0.n_rows;
  arma::uword zMax = z.max();
  //
  arma::mat thetaZ;
  arma::vec sigmaF;
  arma::mat sigmaInt;
  if (logParam) {
    aIntParam2vec(thetaZ, sigmaF, sigmaInt, arma::exp(param), xzDim, zMax);
  } else{
    aIntParam2vec(thetaZ, sigmaF, sigmaInt, param, xzDim, zMax);
  }
  //
  arma::vec y0(n0, fill::zeros);
  arma::vec mse(n0, fill::zeros);
  arma::vec ei(n0, fill::zeros);
  arma::vec ei_1(n0, fill::zeros);
  arma::vec ei_2(n0, fill::zeros);
  aIntNewData(y0, mse, ei, ei_1, ei_2, ei_alpha, min_y, x0, z0, y, x, z, xzDim, 
              mu, invPsi, thetaZ, sigmaF, sigmaInt);
  //
  return List::create(Named("pred") = wrap(y0),
                      Named("mse") = wrap(mse),
                      Named("ei") = wrap(ei),
                      Named("improvement") = wrap(ei_1),
                      Named("uncertainty") = wrap(ei_2)
  );
}

// Independence across Stages

//[[Rcpp::export]]
double lgpNvObjCpp(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget, bool logParam)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu, sigma;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  arma::mat theta;
  //double nugget = 0.;
  if (logParam) {
    nvParam2vec(theta, arma::exp(param), xzDim, zMax);
  } else{
    nvParam2vec(theta, param, xzDim, zMax);  
  }
  nvLogLik(negloglik, psi, invPsi, mu, sigma, nugget, y, x, z, xzDim, theta);
  return negloglik;
}

//[[Rcpp::export]]
Rcpp::List lgpNvModel(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget, bool logParam)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu, sigma;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  arma::mat theta;
  //double nugget = 0.;
  if (logParam) {
    nvParam2vec(theta, arma::exp(param), xzDim, zMax);
  } else{
    nvParam2vec(theta, param, xzDim, zMax);  
  }
  nvLogLik(negloglik, psi, invPsi, mu, sigma, nugget, y, x, z, xzDim, theta);
  //
  return List::create(Named("mu") = wrap(mu),
                      Named("sigma2") = wrap(sigma),
                      Named("theta") = wrap(theta),
                      Named("psi") = wrap(psi),
                      Named("invPsi") = wrap(invPsi),
                      Named("negloglik") = wrap(negloglik),
                      Named("nugget") = wrap(nugget),
                      Named("vecParams") = wrap(param)
  );
}

//[[Rcpp::export]]
Rcpp::List lgpNvPred(arma::mat x0, arma::uvec z0, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim,
                     arma::rowvec param, arma::mat invPsi, double mu, double sigma, double ei_alpha, double min_y, bool logParam)
{
  arma::uword n0 = x0.n_rows;
  arma::uword zMax = z.max();
  //
  arma::mat theta;
  if (logParam) {
    nvParam2vec(theta, arma::exp(param), xzDim, zMax);
  } else{
    nvParam2vec(theta, param, xzDim, zMax);  
  }
  //
  arma::vec y0(n0, fill::zeros);
  arma::vec mse(n0, fill::zeros);
  arma::vec ei(n0, fill::zeros);
  arma::vec ei_1(n0, fill::zeros);
  arma::vec ei_2(n0, fill::zeros);
  nvNewData(y0, mse, ei, ei_1, ei_2, ei_alpha, min_y, x0, z0, y, x, z, xzDim, 
            mu, sigma, invPsi, theta);
  //
  return List::create(Named("pred") = wrap(y0),
                      Named("mse") = wrap(mse),
                      Named("ei") = wrap(ei),
                      Named("improvement") = wrap(ei_1),
                      Named("uncertainty") = wrap(ei_2)
  );
}





