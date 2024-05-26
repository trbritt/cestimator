#include "evaluators.hpp"
std::vector <double> matrix2vector(MatrixXd& x) {
   MatrixXd xx = x;
   int      N, M;
   N = x.rows();
   M = x.cols();

   std::vector <double> data;
   xx.resize(1, N * M);
   for (int i = 0; i < N * M; ++i) {
      data.push_back(xx(i));
   }
   return data;
}

namespace Cestimator {
MatrixXd Evaluator::data, Evaluator::sigma;
VectorXd Evaluator::mu;
int      Evaluator::N, Evaluator::T;

MatrixXd Evaluator::generator(const std::optional <int> nu) {
   double   dof       = static_cast <double>(nu.value_or(1));
   MatrixXd inv_sigma = sigma.llt().solve(MatrixXd::Identity(N, N));
   //create z-score with inv_sigma * (data-mu)
   MatrixXd z_c = data.colwise() - mu;
   switch (static_cast<int>(dist)) {

   case 0:
      return z_c.unaryExpr([&dof](double x) {
                      return exp(-0.5 * x);
                   });
   case 1:
      return z_c.unaryExpr([&dof](double x) {
                      return (x >= 0 && x <= 1) ? gamma(static_cast <double>(N) / 2 + 1) * pow(M_PI, -static_cast <double>(N) / 2)  : 0.0;
                   });

   case 2:
      return z_c.unaryExpr([&dof](double x) {
                      return pow((1 + x / dof), -0.5 * (dof + static_cast <double>(N))) * gamma(0.5 * (dof + static_cast <double>(N))) / (gamma(0.5 * dof) * pow(dof * M_PI, 0.5 * N));
                   });
   }
   ;
};
MatrixXd Evaluator::finite_difference(MatrixXd& gen, int order) {
   //we compute finite differences along ech dimension of the generated array, and evaluate at 0. the fdcoeffF function determines the coefficient of every irregularly spaced data point in the finite differences approximation up to order k, so we just sum the values then
   MatrixXd coeffs = MatrixXd::Zero(N, T);
   for (int n = 0; n < N; ++n) {
      coeffs.row(n) = fdcoeffF(order, 0.0, data.row(n));
   }
   return gen.cwiseProduct(coeffs).rowwise().sum();
}

//taken from https://faculty.washington.edu/rjl/fdmbook/matlab/fdcoeffF.m
MatrixXd Evaluator::fdcoeffF(int k, double xbar, const VectorXd &x, bool fullC ) {
    int n = x.size() - 1;
    if (k > n) {
        throw std::invalid_argument("*** len(x) must be larger than k");
    }

    int m = k;  // for consistency with Fornberg's notation
    double c1 = 1.0;
    double c4 = x(0) - xbar;
    MatrixXd C = MatrixXd::Zero(n + 1, m + 1);
    C(0, 0) = 1.0;

    for (int i = 1; i <= n; ++i) {
        int mn = std::min(i, m);
        double c2 = 1.0;
        double c5 = c4;
        c4 = x(i) - xbar;

        for (int j = 0; j < i; ++j) {
            double c3 = x(i) - x(j);
            c2 = c2 * c3;
            if (j == i - 1) {
                for (int s = mn; s > 0; --s) {
                    C(i, s) = c1 * (s * C(i - 1, s - 1) - c5 * C(i - 1, s)) / c2;
                }
                C(i, 0) = -c1 * c5 * C(i - 1, 0) / c2;
            }
            for (int s = mn; s > 0; --s) {
                C(j, s) = (c4 * C(j, s) - s * C(j, s - 1)) / c3;
            }
            C(j, 0) = c4 * C(j, 0) / c3;
        }
        c1 = c2;
    }

    if (fullC) {
        return C;
    } else {
        return C.col(m); // last column of C
    }
}
}