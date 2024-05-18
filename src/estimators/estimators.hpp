/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: estimators.hpp
    license: gpl_v3

    this is the header for all of the estimators currently implemented in the programme

*/
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>

using namespace Eigen;

typedef std::tuple<Vector2d, MatrixXd> estimator_result;

estimator_result non_parametric_estimation(MatrixXd x);

estimator_result shrinkage_estimation(MatrixXd x);

estimator_result hubert_M(MatrixXd x);