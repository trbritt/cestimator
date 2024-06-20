#include "sde.hpp"

VectorXd Cestimator::SDE::Propagator::Exact::next(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
   return this->model()->exact_stepping(t, dt, x, dZ);
}

VectorXd Cestimator::SDE::Propagator::Euler::next(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
   VectorXd xp = x.array() + this->model()->mu(x, t).array() * dt + this->model()->sigma(x, t).array() * sqrt(dt) * dZ.array();
   return this->model()->positivity() ? xp.array().max(0.0) : xp;
}

VectorXd Cestimator::SDE::Propagator::Milstein::next(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
   VectorXd xp = x.array() + this->model()->mu(x, t).array() * dt +
                 this->model()->sigma(x, t).array() * sqrt(dt) * dZ.array() +
                 0.5 * this->model()->sigma(x, t).array() * this->model()->dsigma_dX(x, t).array() * (dZ.cwiseProduct(dZ).array() - 1) * dt;
   return this->model()->positivity() ? xp.array().max(0.0) : xp;
}
