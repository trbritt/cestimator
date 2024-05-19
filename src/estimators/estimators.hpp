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

namespace Estimator{

    typedef std::tuple<VectorXd, MatrixXd> Result;

    Result non_parametric(MatrixXd x);

    Result shrinkage(MatrixXd x);

    Result maximum_likelihood(MatrixXd x);

    Result robust(MatrixXd x);
}