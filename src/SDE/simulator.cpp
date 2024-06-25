#include "sde.hpp"
#include <random>

MatrixXd generate_normal_matrix(int rows, int cols, double mean = 0.0, double stddev = 1.0) {
   // Random number generator
   std::random_device          rd;
   std::mt19937                gen(rd());
   std::normal_distribution <> d(mean, stddev);

   // Generate matrix
   MatrixXd norms(rows, cols);
   for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
         norms(i, j) = d(gen);
      }
   }

   return norms;
};

namespace Cestimator {
namespace SDE {
MatrixXd Simulator::_initialise_path(const int num_times) {
   MatrixXd p = MatrixXd::Zero(_n_paths, num_times); //to be consistent with other data, each row is a path, and realisations of this path are along second dimension
   p(all, 0).array() = _S0;
   return p;
}

MatrixXd Simulator::_simulate_substep() {
   MatrixXd path   = this->_initialise_path(_M * _substep + 1);
   MatrixXd norms  = generate_normal_matrix(_n_paths, _M * _substep);
   double   dt_sub = _dt / _substep;
   for (int i = 0; i < _M * _substep; ++i) {
      path(all, i + 1) = _pmodel.get()->next(i * dt_sub, dt_sub, path(all, i), norms(all, i));
   }
   return path(all, seq(0, last, _substep));
}

MatrixXd Simulator::simulate_paths() {

   if (_substep > 1 && _pmodel.get()->simulation_method() == "exact") {
      return this->_simulate_substep();
   }
   MatrixXd path  = this->_initialise_path(_M + 1);
   MatrixXd norms = generate_normal_matrix(_n_paths, _M);
   for (int i = 0; i < _M; ++i) {
      path(all, i + 1) = _pmodel.get()->next(i * _dt, _dt, path(all, i), norms(all, i));
   }
   return path(all, seq(0, last, _substep));
}
}
}
