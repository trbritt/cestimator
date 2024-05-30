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
};

namespace Cestimator {
MatrixXd Evaluator::data, Evaluator::sigma;
VectorXd Evaluator::mu;
int      Evaluator::N, Evaluator::T;

VectorXd Evaluator::pad_vector(const VectorXd& d, int pad) {
   int      n = d.size();
   VectorXd padded(n + 2 * pad);
   for (int i = 0; i < pad; ++i) {
      padded(i)           = d(pad - i - 1);
      padded(n + pad + i) = d(n - i - 1);
   }
   for (int i = 0; i < n; ++i) {
      padded(pad + i) = d(i);
   }
   return padded;
};

VectorXd Evaluator::lanczos_kernel(int winlen) {
   int      m     = (winlen - 1) / 2;
   double   denom = m * (m + 1.0) * (2 * m + 1.0);
   VectorXd k     = VectorXd::LinSpaced(m, 1, m);
   VectorXd f     = 3 * k / denom;

   VectorXd differentiator(winlen);
   differentiator.segment(0, m)     = -f.reverse();
   differentiator(m)                = 0.0;
   differentiator.segment(m + 1, m) = f;

   return differentiator;
};

MatrixXd Evaluator::generator(const std::optional <int> nu) {
   double   dof       = static_cast <double>(nu.value_or(1));
   MatrixXd inv_sigma = sigma.llt().solve(MatrixXd::Identity(N, N));
   //create z-score with inv_sigma * (data-mu)
   MatrixXd z_c = data.colwise() - mu;
   switch (static_cast <int>(dist)) {
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
VectorXd Evaluator::differentiate(const MatrixXd& gen, int dimension, int winlen) {
   //we compute finite differences along ech dimension of the generated array, and evaluate at 0. the fdcoeffF function determines the coefficient of every irregularly spaced data point in the finite differences approximation up to order k, so we just sum the values then
   if ((winlen % 2) != 1) {
      throw std::invalid_argument("`winlen` must be an odd number.");
   }
   int      pad           = winlen / 2;
   VectorXd padded_data   = pad_vector(gen.row(dimension), pad);
   VectorXd padded_coords = pad_vector(data.row(dimension), pad);
   VectorXd kernel        = lanczos_kernel(winlen);

   int      n = padded_data.size();
   VectorXd gprime(n - 2 * pad);

   for (int i = pad; i < n - pad; ++i) {
      double sum = 0.0;
      for (int j = -pad; j <= pad; ++j) {
         double dx = padded_coords(i + j) - padded_coords(i);
         sum += padded_data(i + j) * kernel(pad + j) / dx;
      }
      gprime(i - pad) = sum;
   }
   return gprime;
//    int    closest_to_origin = -1;
//    double min_dist          = std::numeric_limits <double>::max();
//    for (int t = 0; t < T; ++t) {
//       double dist = data.col(t).norm();
//       if (dist < min_dist) {
//          min_dist          = dist;
//          closest_to_origin = t;
//       }
//    }
//    VectorXd gprime_at_origin = gprime.col(closest_to_origin);
//    return gprime_at_origin;
}
}
