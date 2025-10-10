// HEADER
void ahgpParam2vec(double &theta, double &sigma, arma::mat &thetaZ, arma::vec &sigmaZ, 
                  const arma::rowvec &param, const arma::uword &xzDim, const arma::uword &zMax);
                  
double ahgpCorrKern(const arma::rowvec &xi, const arma::rowvec &xj, const arma::uword &zi, const arma::uword &zj,
                   const arma::uword &xzDim, const double &theta, const double &sigma, 
                   const arma::mat &thetaZ, const arma::vec &sigmaZ);

void ahgpCorrMat(arma::mat &psi, const arma::mat &x, const arma::uvec &z, const arma::uword &xzDim, 
                const double &theta, const double &sigma, const arma::mat &thetaZ, const arma::vec &sigmaZ);

void ahgpCorrVecs(arma::mat &phi, const arma::mat &x0, const arma::uvec &z0, 
                 const arma::mat &x, const arma::uvec &z, const arma::uword &xzDim, 
                 const double &theta, const double &sigma, const arma::mat &thetaZ, const arma::vec &sigmaZ);

void ahgpCorrVAR(arma::vec &predvar, const arma::uvec &z0, const double &sigma, const arma::vec &sigmaZ);
  
void ahgpLogLik(double &negloglik, arma::mat &psi, arma::mat &invPsi, double &mu, double &nugget,
               const arma::vec &y, const arma::mat &x, const arma::uvec &z, const arma::uword &xzDim, 
               const double &theta, const double &sigma, const arma::mat &thetaZ, const arma::vec &sigmaZ);

void ahgpNewData(arma::vec &y0, arma::vec &mse, arma::vec &ei, arma::vec &ei_1, arma::vec &ei_2, double &ei_alpha, double &min_y,
                const arma::mat &x0, const arma::uvec &z0, const arma::vec &y, const arma::mat &x, const arma::uvec &z, const arma::uword &xzDim, 
                double &mu, arma::mat &invPsi, const double &theta, const double &sigma, const arma::mat &thetaZ, const arma::vec &sigmaZ);

// BODY
void ahgpParam2vec(double &theta, double &sigma, arma::mat &thetaZ, arma::vec &sigmaZ, 
                  const arma::rowvec &param, const arma::uword &xzDim, const arma::uword &zMax)
{
  /* 
   START ASSIGN PARAMETER POSITION 
   */
  thetaZ.set_size(zMax, xzDim); 
  sigmaZ.set_size(zMax);
  /* 
   Parameters for Continuous variables
   */
  theta = param(0);
  arma::uword n_thetaZ = zMax*xzDim;
  thetaZ = arma::reshape(param.subvec(1, n_thetaZ), zMax, xzDim);
  /* 
   Parameters for Variances
   */
  sigma = param(n_thetaZ + 1);
  arma::uword ct = n_thetaZ + 2;
  for (arma::uword u = 0; u < zMax; u++) { sigmaZ(u) = param(ct); ct++; }
}

// CORRELATION KERNEL OF GAUSSIAN PROCESS
double ahgpCorrKern(const arma::rowvec &xi, const arma::rowvec &xj, const arma::uword &zi, const arma::uword &zj,
                   const arma::uword &xzDim, const double &theta, const double &sigma, 
                   const arma::mat &thetaZ, const arma::vec &sigmaZ) 
{
  double zzi = (double)zi;
  double zzj = (double)zj;
  double corrZ = sigma*sigma*std::exp(-(1.0)*theta*(zzi - zzj)*(zzi - zzj));
  /* corr X*/
  arma::rowvec xDiff = xi - xj;
  arma::uword zComm = 0;
  if (zi > zj) { zComm = zj; } else { zComm = zi; }
  double corrX = 0.0;
  for (arma::uword i = 0; i < zComm; i++) {
    arma::rowvec xdtmp = xDiff.subvec(i*xzDim, (i+1)*xzDim - 1);
    corrX += sigmaZ(i)*sigmaZ(i)*std::exp(-(1.0)*arma::accu(thetaZ.row(i) % xdtmp));
  }
  double val = corrZ + corrX;
  return val;
}


void ahgpCorrMat(arma::mat &psi, const arma::mat &x, const arma::uvec &z, const arma::uword &xzDim, 
                const double &theta, const double &sigma, const arma::mat &thetaZ, const arma::vec &sigmaZ)
{
  arma::uword n = x.n_rows;
  for (uword i = 0; i < n; i++) {
    for (uword j = 0; j < i; j++) {
      arma::rowvec xi = x.row(i); 
      arma::rowvec xj = x.row(j);
      arma::uword zi = z(i);
      arma::uword zj = z(j);
      double ker = ahgpCorrKern(xi, xj, zi, zj, xzDim, theta, sigma, thetaZ, sigmaZ);
      psi(i, j) = ker;
      psi(j, i) = ker;
    }
  }
}

void ahgpCorrVecs(arma::mat &phi, const arma::mat &x0, const arma::uvec &z0, 
                 const arma::mat &x, const arma::uvec &z, const arma::uword &xzDim, 
                 const double &theta, const double &sigma, const arma::mat &thetaZ, const arma::vec &sigmaZ)
{
  arma::uword n = x.n_rows;
  arma::uword n0 = x0.n_rows;
  for (uword j = 0; j < n0; j++) {
    arma::rowvec x0j = x0.row(j);
    arma::uword z0j = z0(j);
    for (uword i = 0; i < n; i++) {
      arma::rowvec xi = x.row(i); 
      arma::uword zi = z(i); 
      double ker = ahgpCorrKern(xi, x0j, zi, z0j, xzDim, theta, sigma, thetaZ, sigmaZ);
      phi(i, j) = ker;
    }
  }
}

void ahgpCorrVAR(arma::vec &predvar, const arma::uvec &z0, const double &sigma, const arma::vec &sigmaZ)
{
  arma::uword n0 = z0.n_elem;
  for (uword j = 0; j < n0; j++) {
    predvar(j) = sigma;
    arma::uword z0j = z0(j);
    for (arma::uword i = 0; i < z0j; i++) {
      predvar(j) += sigmaZ(i);
    }
  }
}

void ahgpLogLik(double &negloglik, arma::mat &psi, arma::mat &invPsi, double &mu, double &nugget,
               const arma::vec &y, const arma::mat &x, const arma::uvec &z, const arma::uword &xzDim, 
               const double &theta, const double &sigma, const arma::mat &thetaZ, const arma::vec &sigmaZ)
{
  arma::uword n = y.n_elem;
  arma::vec onevec(n, fill::ones);
  ahgpCorrMat(psi, x, z, xzDim, theta, sigma, thetaZ, sigmaZ);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, psi);
  arma::mat eyemat(n, n, fill::eye);
  double checkCond = std::abs(eigval.max()) - 1e8*std::abs(eigval.min());
  if ((nugget == 0) & (checkCond >= 0)) {
    nugget = checkCond/(1e8 - 1);
  }
  psi += nugget*eyemat;
  double detPsi;
  double signDetPsi;
  bool invSucc;
  invSucc = arma::inv_sympd(invPsi, psi);
  arma::log_det(detPsi, signDetPsi, psi);
  //if (std::isfinite(detPsi) & (signDetPsi >= 0)) 
  if (invSucc) {
    double yPsiY = arma::as_scalar(y.t()*invPsi*y);
    double onePsiY = arma::as_scalar(onevec.t()*invPsi*y);
    double onePsiOne = arma::as_scalar(onevec.t()*invPsi*onevec);
    mu = onePsiY/onePsiOne;
    negloglik = (-1.0)*(-0.5)*(detPsi + yPsiY - (onePsiY*onePsiY)/onePsiOne);
  } else {
    negloglik = 1e20;
  }
}


void ahgpNewData(arma::vec &y0, arma::vec &mse, arma::vec &ei, arma::vec &ei_1, arma::vec &ei_2, double &ei_alpha, double &min_y,
                const arma::mat &x0, const arma::uvec &z0, const arma::vec &y, const arma::mat &x, const arma::uvec &z, const arma::uword &xzDim, 
                double &mu, arma::mat &invPsi, const double &theta, const double &sigma, const arma::mat &thetaZ, const arma::vec &sigmaZ)
{
  arma::uword n = x.n_rows;
  arma::uword n0 = x0.n_rows;
  arma::mat phi(n, n0, fill::zeros);
  ahgpCorrVecs(phi, x0, z0, x, z, xzDim, theta, sigma, thetaZ, sigmaZ);
  arma::vec predvar(n0, fill::zeros);
  ahgpCorrVAR(predvar, z0, sigma, sigmaZ);
  arma::vec onevec(n, fill::ones);
  arma::vec resid = y - mu*onevec; 
  arma::vec psiinvresid = invPsi*resid;
  for (uword j = 0; j < n0; j++) {
    y0(j) = mu + arma::as_scalar(phi.col(j).t()*psiinvresid);
    mse(j) = std::abs(predvar(j) - arma::as_scalar(phi.col(j).t()*invPsi*phi.col(j))) + datum::eps;
  }
  // Compute expected improvement
  //double min_val = arma::min(y);
  arma::vec rmse = arma::sqrt(mse);
  arma::vec yd = min_y - y0;
  // The improvement part
  ei_1 = yd % (.5 + .5*arma::erf((1./std::sqrt(2.))*(yd/rmse)));
  // The uncertainty part
  ei_2 = (rmse/std::sqrt(2.*datum::pi)) % arma::exp(-.5*(yd % yd)/mse);
  // The EI value
  ei = 2.*(ei_alpha*ei_1 + (1. - ei_alpha)*ei_2);
  ei.elem( arma::find(ei <= .0) ).fill(datum::eps);  
}

