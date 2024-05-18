/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: nonparametric.cpp
    license: gpl_v3

    this is the main code for the nonparametric estimator

*/

#include "estimators.hpp"
#include "../utils/utils.hpp"

estimator_result non_parametric_estimation(MatrixXd x){
    Vector2d mu = mean(x);
    
    MatrixXd sigma = covariance(x);

    return std::make_tuple(mu, sigma);
}