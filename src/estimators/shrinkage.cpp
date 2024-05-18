#include "estimators.hpp"
#include "../utils/utils.hpp"

estimator_result shrinkage_estimation(MatrixXd x){

    const int N = x.rows();
    const int T = x.cols();

    Vector2d mu_no_par = mean(x);
    
    MatrixXd sigma_no_par = covariance(x);

    MatrixXd b = MatrixXd::Zero(N, 1);

    EigenSolver<MatrixXd> sigma_solver(sigma_no_par);

    VectorXd lambda_hat = sigma_solver.eigenvalues().real();

    double a = (1/T) * (lambda_hat.sum() - 2*lambda_hat.maxCoeff()) / (mu_no_par.transpose() * mu_no_par);
    a = std::fmax(0.0, std::fmin(a, 1.0));

    VectorXd mu_shr = (1-a) * mu_no_par + a * b;

    MatrixXd C = lambda_hat.mean() * MatrixXd::Identity(N,N);

    // now we compute the optimal weight
    double num = 0;
    MatrixXd tmp;
    for (int t = 0; t< T; ++t){
        tmp = x(all, t) * x(all,t).transpose() - sigma_no_par;
        tmp *= tmp;
        num += (1/T) * tmp.trace(); 
    };
    
    tmp = (sigma_no_par - C) * (sigma_no_par - C);
    double denom = tmp.trace();

    a = (1/T) * num / denom;
    a = std::fmax(0.0, std::fmin(a, 1.0));

    MatrixXd sigma_shr = (1-a) * sigma_no_par + a*C;
    return std::make_tuple(mu_shr, sigma_shr);
}