#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>


typedef std::tuple<Eigen::Vector2d, Eigen::MatrixXd> estimator_result;

estimator_result non_parametric_estimation(Eigen::MatrixXd x);

estimator_result shrinkage_estimation(Eigen::MatrixXd x);

Eigen::MatrixXd covariance(Eigen::MatrixXd x);

Eigen::Vector2d mean(Eigen::MatrixXd x);

Eigen::VectorXd outlier_cutoff(Eigen::VectorXd d, double d0);

estimator_result hubert_M(Eigen::MatrixXd x);