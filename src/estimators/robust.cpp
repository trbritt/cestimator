/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: robust.cpp
    license: gpl_v3

    this is the main code for a robust estimator

*/
#include "estimators.hpp"

Estimator::Result Estimator::robust(MatrixXd x){
    const double tolerance = 1e-6;
    const double error = 1e6;

    const int N = x.rows();
    const int T = x.cols();

    VectorXd w = VectorXd::Ones(T);
    VectorXd Zeros = VectorXd::Zero(N);

    VectorXd mu = Zeros;

    MatrixXd sigma = MatrixXd::Zero(N,N);

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
            mu += w(t)*w(t) * x(all, t).transpose();
        }
        mu /= w.sum();

        sigma = MatrixXd::Zero(N,N);
        for (size_t t=0; t<T; ++t){
            
            tmp = (x(all,t).transpose() - mu);
            sigma += w(t) * w(t) * tmp * tmp.transpose();
        }
        sigma /= (w.transpose() * w);

        inv_sigma = sigma.inverse();
        std::cout << inv_sigma << std::endl;
        break;
        // for (size_t t=0; t<T; ++t){
        //     MatrixXd tmp2 = x(all,t).transpose() - mu;

        //     auto i = tmp2 * inv_sigma * tmp2.transpose();
        //     std::cout << i << std::endl;
        // }

    }

    return std::make_tuple(mu, sigma);
}