/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: utils.hpp
    license: gpl_v3

    this is the header file for all utilities of the project,
    namely anything not theoretically driving the main results,
    but necessary to get there

*/

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>

using namespace Eigen;

MatrixXd covariance(MatrixXd x);

Vector2d mean(MatrixXd x);

VectorXd outlier_cutoff(VectorXd d, double d0);