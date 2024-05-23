#include "estimators.hpp"

namespace Cestimator {
    MatrixXd Estimator::data,Estimator::sigma_no_par;
    VectorXd Estimator::mu_no_par;
    int Estimator::N, Estimator::T;
}