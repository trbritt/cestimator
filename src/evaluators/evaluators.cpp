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
   std::function <double(double)> _generator = [&dof](double x) {
                                                  return 0.0;
                                               };
   switch (dist) {
   case Uniform:
      _generator = [&dof](double x) {
                      return (x >= 0 && x <= 1) ? gamma(static_cast <double>(N) / 2 + 1) * pow(M_PI, -static_cast <double>(N) / 2)  : 0.0;
                   };

   case Normal:
      _generator = [&dof](double x) {
                      return exp(-0.5 * x);
                   };

   case StudentT:
      _generator = [&dof](double x) {
                      return pow((1 + x / dof), -0.5 * (dof + static_cast <double>(N))) * gamma(0.5 * (dof + static_cast <double>(N))) / (gamma(0.5 * dof) * pow(dof * M_PI, 0.5 * N));
                   };

   default:
      _generator = [](double x) {
                      return 1.0;
                   };
   }
   ;
   return z_c.unaryExpr(_generator);
};
MatrixXd Evaluator::finite_difference(MatrixXd& gen, int order) {
   //we compute finite differences along ech dimension of the generated array, and evaluate at 0. the fdcoeffF function determines the coefficient of every irregularly spaced data point in the finite differences approximation up to order k, so we just sum the values then
   MatrixXd coeffs = MatrixXd::Zero(N, T);
   for (int n = 0; n < N; ++n) {
      std::vector <double> flattened = matrix2vector(data);
      std::vector <double> coef      = fdcoeffF(order, 0.0, flattened);
      coeffs.row(n) = MatrixXd::Map(coef.data(), T, 1);
   }
   return gen.cwiseProduct(coeffs).rowwise().sum();
}

//taken from https://faculty.washington.edu/rjl/fdmbook/matlab/fdcoeffF.m
std::vector <double> Evaluator::fdcoeffF(int k, double xbar, const std::vector <double>& x) {
   int n = x.size();
   if (k >= n) {
      throw std::invalid_argument("length(x) must be larger than k");
   }

   int m = k;
   std::vector <std::vector <double> > C(n, std::vector <double>(m + 1, 0.0));
   double c1 = 1.0;
   double c4 = x[0] - xbar;
   C[0][0] = 1.0;

   for (int i = 0; i < n - 1; ++i) {
      int    i1 = i + 1;
      int    mn = std::min(i, m);
      double c2 = 1.0;
      double c5 = c4;
      c4 = x[i1] - xbar;
      for (int j = 0; j <= i; ++j) {
         int    j1 = j + 1;
         double c3 = x[i1] - x[j];
         c2 *= c3;
         if (j == i) {
            for (int s = mn; s >= 1; --s) {
               int s1 = s + 1;
               C[i1][s1] = c1 * (s * C[i1 - 1][s1 - 1] - c5 * C[i1 - 1][s1]) / c2;
            }
            C[i1][0] = -c1 * c5 * C[i1 - 1][0] / c2;
         }
         for (int s = mn; s >= 1; --s) {
            int s1 = s + 1;
            C[j][s1] = (c4 * C[j][s1] - s * C[j][s1 - 1]) / c3;
         }
         C[j][0] = c4 * C[j][0] / c3;
      }
      c1 = c2;
   }

   std::vector <double> c(n);
   for (int i = 0; i < n; ++i) {
      c[i] = C[i][m];
   }

   return c;
}
}
