#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>

using namespace Eigen;

MatrixXd covariance(MatrixXd x);

Vector2d mean(MatrixXd x);

VectorXd outlier_cutoff(VectorXd d, double d0);