/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: shrinkage.cpp
    license: gpl_v3

    this is the main code for the James-Stein shrinkage estimator

*/
#include "estimators.hpp"

int Cestimator::shrinkage::run() noexcept {

    MatrixXd b = MatrixXd::Zero(N, 1);

    EigenSolver<MatrixXd> sigma_solver(sigma_no_par);

    VectorXd lambda_hat = sigma_solver.eigenvalues().real();

    double lambda_hat_sum = lambda_hat.sum();
    double lambda_hat_max = lambda_hat.maxCoeff();
    double mu_no_par_sq2 = mu_no_par.squaredNorm();

    double a = (1.0 / T) * (lambda_hat_sum - 2 * lambda_hat_max) / mu_no_par_sq2;
    a = std::max(0.0, std::min(a, 1.0));

    mu = (1-a) * mu_no_par + a * b;

    MatrixXd C = lambda_hat.mean() * MatrixXd::Identity(N,N);

    // now we compute the optimal weight
    MatrixXd Xcentered;
    MatrixXd sigmacentered = sigma_no_par - C;
    double num = 0;
    for (int t = 0; t< T; ++t){
        Xcentered = data.col(t) * data.col(t).transpose() - sigma_no_par;
        num += (Xcentered.array() * Xcentered.array()).sum();
    };
    
    double denom = (sigmacentered.array() * sigmacentered.array()).sum();

    a = (1/T) * num / denom;
    a = std::fmax(0.0, std::fmin(a, 1.0));

    sigma = (1-a) * sigma_no_par + a*C;
    return 0;
}