#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>

std::tuple<Eigen::Vector2d, Eigen::MatrixXd> non_parametric_estimation(Eigen::MatrixXd x);

std::tuple<Eigen::Vector2d, Eigen::MatrixXd> shrinkage_estimation(Eigen::MatrixXd x);

Eigen::MatrixXd covariance(Eigen::MatrixXd x);

Eigen::Vector2d mean(Eigen::MatrixXd x);