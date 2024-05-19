/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: nonparametric.cpp
    license: gpl_v3

    this is the main code for the nonparametric estimator

*/

#include "estimators.hpp"
#include "../utils/utils.hpp"

Estimator::Result Estimator::non_parametric(MatrixXd x){
    VectorXd mu = Utils::mean(x);
    
    MatrixXd sigma = Utils::covariance(x);

    return std::make_tuple(mu, sigma);
}