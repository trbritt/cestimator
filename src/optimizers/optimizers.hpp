#include <optional>
#include <Eigen/Dense>
#include <iterator>
#include <iostream>
#include <set>
#include <map>
#include <algorithm>

using namespace Eigen;

namespace Cestimator {
class Optimizer {
    public: 
    enum Classifier { Random = 0, Qp = 1};
    Classifier classifier;

    Optimizer(enum Classifier c){
        classifier = c;
    }
    MatrixXd generate_frontier(VectorXd returns, MatrixXd covariance, int n_portfolios);
};
}
