// DECLARE FUNCTIONS
void matrixPrintf(const mat &m);
void vecPrintf(const vec &v);
void rvecPrintf(const rowvec &v);
void autoNugget(arma::mat &psi, double &nugget);
void calcMatrixInv(bool &invSucc, arma::mat &invPsi, arma::mat &psi, double &nugget);

// BODY
void matrixPrintf(const mat &m)
{
  for (uword i = 0; i < m.n_rows; i++) {
    for (uword j = 0; j < m.n_cols; j++) Rprintf("%4.4f\t", m(i,j));
    Rprintf("\n");
  }
	Rprintf("\n\n");
}

void rvecPrintf(const rowvec &v)
{
  for (uword i = 0; i < v.n_elem; i++) Rprintf("%4.4f\t", v(i));
  Rprintf("\n\n");
}

void vecPrintf(const vec &v)
{
  for (uword i = 0; i < v.n_elem; i++) Rprintf("%4.4f\n", v(i));
  Rprintf("\n\n");
}


void autoNugget(arma::mat &psi, double &nugget) 
{
  arma::uword n = psi.n_rows;
  arma::mat eyemat(n, n, fill::eye);
  bool ISSYMPD = psi.is_sympd();
  if (!ISSYMPD && (nugget == 0)) {
    for (arma::uword ng = 0; ng < 101; ng++) {
      nugget = std::exp((100. - (double)(ng))*(-52.)/100.);
      psi += nugget*eyemat;    
      ISSYMPD = psi.is_sympd();
      if (ISSYMPD) { 
        break; 
      }
    }
  }
}

void calcMatrixInv(bool &invSucc, arma::mat &invPsi, arma::mat &psi, double &nugget) 
{
  arma::uword n = psi.n_rows;
  if (nugget == 0) {
    autoNugget(psi, nugget);
  } else {
    arma::mat eyemat(n, n, fill::eye);
    psi += nugget*eyemat;
  }
  /*
  // LU Decomposition Approach
  arma::mat Lmat(n, n, fill::zeros);
  arma::mat invLmat(n, n, fill::zeros);
  bool cholOK = arma::chol(Lmat, psi, "lower");
  Rprintf("%d\n", cholOK);
  if (cholOK) {
    invSucc = arma::inv(invLmat, Lmat);
    invPsi = invLmat.t()*invLmat;
  }
  */
  // Symmetric Positively Definite Approach
  invSucc = arma::inv_sympd(invPsi, psi);
  //Rprintf("%d\n", invSucc);
}

