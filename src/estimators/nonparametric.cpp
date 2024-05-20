/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: nonparametric.cpp
    license: gpl_v3

    this is the main code for the nonparametric estimator

*/

#include "estimators.hpp"
#include "../utils/utils.hpp"

Cestimator::Result Cestimator::non_parametric(MatrixXd x){
    VectorXd mu = Cestimator::Utils::mean(x);
    
    MatrixXd sigma = Cestimator::Utils::covariance(x);

    return std::make_tuple(mu, sigma);
}