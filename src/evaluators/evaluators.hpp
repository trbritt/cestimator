#include <optional>
#include <Eigen/Dense>
#include <iterator>
#include <iostream>

using namespace Eigen;

namespace Cestimator {
class Evaluator {
public:
   enum Distribution { Normal = 0, Uniform = 1, StudentT = 2 };
   Distribution dist;

   static MatrixXd data, sigma;
   static VectorXd mu;
   static int N, T;

   Evaluator(enum Distribution d) {
      dist = d;
   };
   void set_data(MatrixXd& x, VectorXd& m, MatrixXd& s) {
      data  = x;
      mu    = m;
      sigma = s;
      N     = x.rows();
      T     = x.cols();
   };
   MatrixXd generator(const std::optional <int> nu = std::nullopt);
   MatrixXd finite_difference(MatrixXd& gen, int order);
   std::vector <double> fdcoeffF(int k, double xbar, const std::vector <double>& x);
};
}
