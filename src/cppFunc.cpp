
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
double lgpOBObjCpp(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget)
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
  lgpOBParam2vec(tau, theta, param, xzDim, zMax);
  lgpOBLogLik(negloglik, psi, invPsi, mu, sigma2, nugget, y, x, z, xzDim, 
              tau, theta);
  return negloglik;
}

//[[Rcpp::export]]
Rcpp::List lgpOBModel(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget)
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
  lgpOBParam2vec(tau, theta, param, xzDim, zMax);
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
Rcpp::List lgpOBPred(arma::mat x0, arma::umat z0, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim,
                     arma::rowvec param, arma::mat invPsi, double mu, double sigma2, double ei_alpha, double min_y)
{
  arma::uword n0 = x0.n_rows;
  arma::uword zMax = z.max();
  //
  double tau;
  arma::mat theta;
  lgpOBParam2vec(tau, theta, param, xzDim, zMax);
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
double lgpOBNbjCpp(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu, sigma2;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  arma::mat tau, theta;
  //double nugget = 0.;
  lgpNBParam2vec(tau, theta, param, xzDim, zMax);
  lgpNBLogLik(negloglik, psi, invPsi, mu, sigma2, nugget, y, x, z, xzDim, 
              tau, theta);
  return negloglik;
}

//[[Rcpp::export]]
Rcpp::List lgpNBModel(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu, sigma2;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  arma::mat tau, theta;
  //double nugget = 0.;
  lgpNBParam2vec(tau, theta, param, xzDim, zMax);
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
Rcpp::List lgpNBPred(arma::mat x0, arma::umat z0, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim,
                     arma::rowvec param, arma::mat invPsi, double mu, double sigma2, double ei_alpha, double min_y)
{
  arma::uword n0 = x0.n_rows;
  arma::uword zMax = z.max();
  //
  arma::mat tau, theta;
  lgpNBParam2vec(tau, theta, param, xzDim, zMax);
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







/*
// AddHGP
//[[Rcpp::export]]
double ahgpObjCpp(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  double theta, sigma;
  arma::mat thetaZ;
  arma::vec sigmaZ;
  //double nugget = 0.;
  ahgpParam2vec(theta, sigma, thetaZ, sigmaZ, param, xzDim, zMax);
  ahgpLogLik(negloglik, psi, invPsi, mu, nugget, y, x, z, xzDim, 
             theta, sigma, thetaZ, sigmaZ);
  return negloglik;
}

//[[Rcpp::export]]
Rcpp::List ahgpModel(arma::rowvec param, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim, double nugget)
{
  arma::uword n = x.n_rows;
  arma::uword zMax = z.max();
  double negloglik, mu;
  arma::mat psi(n, n, fill::eye); 
  arma::mat invPsi(n, n, fill::eye);
  //
  double theta, sigma;
  arma::mat thetaZ;
  arma::vec sigmaZ;
  //double nugget = 0.;
  ahgpParam2vec(theta, sigma, thetaZ, sigmaZ, param, xzDim, zMax);
  ahgpLogLik(negloglik, psi, invPsi, mu, nugget, y, x, z, xzDim, 
             theta, sigma, thetaZ, sigmaZ);
  //
  return List::create(Named("mu") = wrap(mu),
                      Named("theta") = wrap(theta),
                      Named("thetaZ") = wrap(thetaZ),
                      Named("sigma") = wrap(sigma),
                      Named("sigmaZ") = wrap(sigmaZ),
                      Named("psi") = wrap(psi),
                      Named("invPsi") = wrap(invPsi),
                      Named("negloglik") = wrap(negloglik),
                      Named("nugget") = wrap(nugget),
                      Named("vecParams") = wrap(param)
  );
}

//[[Rcpp::export]]
Rcpp::List ahgpPred(arma::mat x0, arma::umat z0, arma::vec y, arma::mat x, arma::uvec z, arma::uword xzDim,
                    arma::rowvec param, arma::mat invPsi, double mu, double ei_alpha, double min_y)
{
  arma::uword n0 = x0.n_rows;
  arma::uword zMax = z.max();
  //
  double theta, sigma;
  arma::mat thetaZ;
  arma::vec sigmaZ;
  ahgpParam2vec(theta, sigma, thetaZ, sigmaZ, param, xzDim, zMax);
  //
  arma::vec y0(n0, fill::zeros);
  arma::vec mse(n0, fill::zeros);
  arma::vec ei(n0, fill::zeros);
  arma::vec ei_1(n0, fill::zeros);
  arma::vec ei_2(n0, fill::zeros);
  ahgpNewData(y0, mse, ei, ei_1, ei_2, ei_alpha, min_y, x0, z0, y, x, z, xzDim, 
              mu, invPsi, theta, sigma, thetaZ, sigmaZ);
  //
  return List::create(Named("pred") = wrap(y0),
                      Named("mse") = wrap(mse),
                      Named("ei") = wrap(ei),
                      Named("improvement") = wrap(ei_1),
                      Named("uncertainty") = wrap(ei_2)
  );
}
 */




