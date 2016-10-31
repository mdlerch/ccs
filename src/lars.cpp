#include <RcppArmadillo.h>
#include <stdio.h>


using namespace arma;

// [[Rcpp::export]]
arma::mat larscpp(arma::mat x, arma::mat beta, arma::vec y, int maxk) {

    int p = x.n_cols;
    int n = x.n_rows;
    int j;

    arma::vec cvec(p);
    arma::vec cveca(p);
    arma::vec mu(n);
    arma::uvec active(n);
    arma::mat XA(1, 1);

    int i = 0;
    uword nv = 0;
    while (i < maxk & nv < p)
    {
        ++i;

        // 1. find next variable to add
        cvec = x.t() * (y - mu);
        cveca = abs(cvec);

        double cmax = 0;
        for (int k = 0; k < p; ++k)
        {
            if (cveca(k) > cmax)
            {
                j = k;
                cmax = cveca(k);
            }
        }

        // 2. find unit-vector of equal projection
        active(nv) = j;
        nv++;
        vec signs = sign(cvec);
        uvec activeshort = active.subvec(0, nv);
        arma::mat XA = x.cols(activeshort);

    }

    return XA;
}
