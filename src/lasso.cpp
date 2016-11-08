#include <RcppArmadillo.h>
#include <stdio.h>

using namespace arma;

// [[Rcpp::export]]
arma::mat lassocpp(arma::mat x, arma::mat beta, arma::vec y, int maxk) {

    int p = x.n_cols;
    int n = x.n_rows;
    int j, out;
    bool lass = false;

    double gam;

    vec cvec(p);
    vec cveca(p);
    ivec cvecs(p);
    uvec active(p);
    uvec inactive(p);
    mat XA(n, p);

    vec mu = x * beta.row(0).t();
    int idx;

    for (int i = 0; i < p; ++i)
    {
        inactive(i) = i;
    }


    int i = 0;
    uword nv = 0;
    while (i < maxk & nv < p)
    {
        ++i;

        // 1. find next variable to add
        cvec = x.t() * (y - mu);
        cveca = abs(cvec);

        double cmax = 0;
        if (lass)
        {
            for (int k = 0; k < p; ++k)
            {
                if (k != out)
                {
                    if (cveca(k) > cmax)
                    {
                        cmax = cveca(k);
                    }
                    if (cvec(k) > 0)
                    {
                        cvecs(k) = 1;
                    } else {
                        cvecs(k) = -1;
                    }
                }
            }
            uvec tempactive(p);
            int tempidx = 0;
            for (int k = 0; k <= nv; ++k)
            {
                if (active(k) != out)
                {
                    tempactive(tempidx) = active(k);
                    tempidx++;
                }
            }
            active = tempactive;
        } else {
            for (int k = 0; k < (p - nv); ++k)
            {
                idx = inactive(k);
                if (cveca(idx) > cmax)
                {
                    j = idx;
                    cmax = cveca(idx);
                }
                if (cvec(idx) > 0)
                {
                    cvecs(idx) = 1;
                } else {
                    cvecs(idx) = -1;
                }
            }
            cout << "in: " << j << endl;
            active(nv) = j;
            nv++;
        }


        // build inactive
        idx = 0;
        bool skip;
        for (int k = 0; k < p; ++k)
        {
            skip = false;
            for (int h = 0; h < nv; ++h)
            {
                if (k == active(h))
                {
                    skip = true;
                    break;
                }
            }
            if (!skip)
            {
                inactive(idx) = k;
                idx++;
            }
        }


        // 2. find unit-vector of equal projection
        for (int k = 0; k < nv; ++k)
        {
            XA.col(k) = x.col(active(k)) * cvecs(active(k));
        }
        mat gA = XA.cols(0, nv - 1).t() * XA.cols(0, nv - 1);
        vec one(nv, fill::ones);
        mat AA = 1 / sqrt(one.t() * pinv(gA) * one);
        vec w = AA(0, 0) * (pinv(gA) * one);
        vec u = XA.cols(0, nv - 1) * w;

        // 3. increment
        vec a = x.t() * u;
        gam = cmax / AA(0, 0);
        for (int k = 0; k < (p - nv); ++k)
        {
            idx = inactive(k);
            double temp = (cmax - cvec(idx)) / (AA(0, 0) - a(idx));
            if (temp > 0 & temp < gam)
            {
                gam = temp;
            }
            temp = (cmax + cvec(idx)) / (AA(0, 0) + a(idx));
            if (temp > 0 & temp < gam)
            {
                gam = temp;
            }
        }
        double gamj;
        lass = false;
        for (int k = 0; k < nv; ++k)
        {
            gamj = - beta(i - 1, active(k)) / (w(k) * cvecs(active(k)));
            if (gamj > 0 & gamj < gam)
            {
                gam = gamj;
                out = active(k);
                lass = true;
                cout << "lasso" << endl; 
            }
        }

        cout << "gamma: " << gam << endl;;


        for (int k = 0; k < nv; ++k)
        {
            beta(i, active(k)) = beta(i - 1, active(k)) + gam * w(k) * cvecs(active(k));
        }
        mu = mu + gam * u;

        if (lass)
        {
            nv--;
        }

        // cout << beta.row(i) << endl;
    }

    return beta.rows(0, i);
}
