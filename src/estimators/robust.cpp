/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: robust.cpp
    license: gpl_v3

    this is the main code for a robust estimator

*/
#include "estimators.hpp"

int Cestimator::robust::run() noexcept {
    const double tolerance = 1e-6;
    double error = 1e6;

    VectorXd w = VectorXd::Ones(T);
    VectorXd Zeros = VectorXd::Zero(N);

    mu = Zeros;

    sigma = MatrixXd::Zero(N,N);

    double d0 = sqrt(N) + sqrt(2);
    VectorXd mu_old;
    MatrixXd sigma_old;
    MatrixXd inv_sigma;
    MatrixXd tmp;
    VectorXd d = VectorXd::Ones(T);

    while (error > tolerance){
        mu_old  = mu;
        sigma_old = sigma;

        mu = Zeros;
        for (size_t t=0; t<T; ++t){
            mu += w(t)*w(t) * data(all, t).transpose();
        }
        mu /= w.sum();

        sigma = MatrixXd::Zero(N,N);
        for (size_t t=0; t<T; ++t){
            
            tmp = (data(all,t) - mu);
            sigma += w(t) * w(t) * tmp * tmp.transpose();
        }
        sigma /= (w.transpose() * w);

        inv_sigma = sigma.inverse();
        for (size_t t=0; t<T; ++t){
            MatrixXd tmp2 = data(all,t) - mu;
            //below is an annoying typecast issue, the matrix products work out to a scalar value, but it is still of the type MatrixXd, so to actually get the value we take the (1,1) element explicitly, despite it being the only element            
            d(t) = sqrt((tmp2.transpose() * inv_sigma * tmp2)(0));
        }
        w = Cestimator::Utils::outlier_cutoff(d, d0);
        error = ((sigma-sigma_old)*(sigma-sigma_old) + (mu-mu_old) * (mu-mu_old).transpose()).trace();
    }

    return 0;
}