/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: utils.cpp
    license: gpl_v3

    this is the implementation for all utilities of the project,
    namely anything not theoretically driving the main results,
    but necessary to get there

*/

#include "utils.hpp"

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