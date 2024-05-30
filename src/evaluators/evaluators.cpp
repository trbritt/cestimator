#include "evaluators.hpp"

bool comparator(const std::pair <double, int>& a, const std::pair <double, int>& b) {
   return a.first < b.first;
}

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
   // we only consider those points immediately around the area of differentiation, then get the unique values, and sort them to stabilze convolution
   std::vector <double> dists(T, std::numeric_limits <double>::max());
   for (int t = 0; t < T; ++t) {
      dists[t] = std::sqrt(std::pow(data.col(t)(dimension), 2));
   }
   std::map <double, int> element2index;
   for (int t = 0; t < T; ++t) {
      if (element2index.find(dists[t]) == element2index.end()) {
         element2index[dists[t]] = t;
      }
   }
   std::vector <std::pair <double, int> > sorted_dists(element2index.begin(), element2index.end());
   std::sort(sorted_dists.begin(), sorted_dists.end(), comparator);   //default behaviour is to use .first. but we are explicit just in case

   //this now puts the points sorted by distance, we take the top points closest to the origin

   int n_neighbours = std::min(
      static_cast <int>(element2index.size()),
      std::min(T, static_cast <int>(10 * winlen))
      );
   std::vector <std::pair <double, int> > sorted_coords(n_neighbours);

   for (int i = 0; i < n_neighbours; ++i) {
      sorted_coords[i] = std::make_pair(data(dimension, sorted_dists[i].second), sorted_dists[i].second);
   }
   std::sort(sorted_coords.begin(), sorted_coords.end(), comparator);

   // but now we want to sort the original coordinates of the data at these indices
   VectorXd u_s_data(n_neighbours);
   VectorXi u_s_indices(n_neighbours);
   double   min_val = std::numeric_limits <double>::max();
   int      min_idx = -1;
   for (int i = 0; i < n_neighbours; ++i) {
      u_s_data(i)    = sorted_coords[i].first;
      u_s_indices(i) = sorted_coords[i].second;
      if (u_s_data(i) < min_idx) {
         min_val = u_s_data(i);
         min_idx = i;
      }
   }
   int      pad           = winlen / 2;
   VectorXd padded_data   = pad_vector(gen.row(dimension)(u_s_indices), pad);
   VectorXd padded_coords = pad_vector(u_s_data, pad);
   VectorXd kernel        = lanczos_kernel(winlen);

   int      n = padded_data.size();
   VectorXd gprime(n - 2 * pad);

   for (int i = pad; i < n - pad; ++i) {
      double sum = 0.0;
      for (int j = -pad; j <= pad; ++j) {
         double dx = padded_coords(i + j) - padded_coords(i);
         if (dx < 1e-6) {
            continue;
         }
         sum += padded_data(i + j) * kernel(pad + j) / dx;
      }
      gprime(i - pad) = sum;
   }
   return gprime;
}
}
