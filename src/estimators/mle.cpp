/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: mle.cpp
    license: gpl_v3

    this is the main code for the maximum likelihood estimator

*/
#include <vector>
#include <iterator>
#include "estimators.hpp"
#include "../utils/utils.hpp"

std::vector<double> matrix2vector(MatrixXd x){
    MatrixXd xx = x;
    int N, M;
    N = x.rows();
    M = x.cols();

    std::vector<double> data;
    xx.resize(1, N*M);
    for(int i=0; i<N*M; ++i){
        data.push_back(xx(i));
    }
    return data;
}

Estimator::Result Estimator::maximum_likelihood(MatrixXd x){

    Vector2d mu = Utils::mean(x);
    
    MatrixXd sigma = Utils::covariance(x);

    std::vector<double> quant = Utils::quantile(matrix2vector(x), {0.75, 0.25});

    const double tolerance = 0.01*(quant[1]-quant[0]);
    const int nus[6] = {1, 2, 4, 7, 12, 20};

    
    return std::make_tuple(mu, sigma);
}