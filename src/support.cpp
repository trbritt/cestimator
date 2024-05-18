#include "support.hpp"

using namespace Eigen;

Vector2d mean(MatrixXd x){
    return x.rowwise().mean();
}

MatrixXd covariance(MatrixXd x){
    MatrixXd tmp = x.transpose();
    MatrixXd centered = tmp.rowwise() - tmp.colwise().mean();
    return (tmp.adjoint() * centered) / (double(tmp.rows() - 1));
}

VectorXd outlier_cutoff(VectorXd d, double d0){
    Matrix<int, -1, 1> index{d.size()};
    for (size_t i=0; i < d.size(); ++i){
        index(i) = d(i) <= d0 ? 1 : 0;
    }
    const double b_2 = 1.25;
    VectorXd omega{d.size()};
    for (size_t i = 0; i < d.size(); ++i){
        if (index(i)) {
            omega(i) = 1;
        } else {
            omega(i) = -0.5*((d(i)-d0)*(d(i)-d0) / (b_2*b_2));
        }
    }
    return omega;
}

estimator_result non_parametric_estimation(MatrixXd x){
    Vector2d mu = mean(x);
    
    MatrixXd sigma = covariance(x);

    return std::make_tuple(mu, sigma);
}

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

estimator_result hubert_M(MatrixXd x){
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