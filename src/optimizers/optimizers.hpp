#include <optional>
#include <Eigen/Dense>
#include <iterator>
#include <iostream>

using namespace Eigen;

namespace Cestimator {
class Optimizer {
public:
   enum Classifier { Random = 0, Qp = 1 };
   Classifier classifier;
   std::vector <std::pair <double, double> > profile;

   Optimizer(enum Classifier c) {
      classifier = c;
   }

   MatrixXd generate_frontier(VectorXd returns, MatrixXd covariance, int n_portfolios);
};
}
