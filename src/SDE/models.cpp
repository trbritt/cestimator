#include "sde.hpp"

VectorXd Cestimator::SDE::Model::mu(const VectorXd& x, double t) {
   switch (_type) {
   case (Brownian):
      return _params(0) * (x.array() > -10000).cast <double>().matrix();
   }
}

VectorXd Cestimator::SDE::Model::sigma(const VectorXd& x, double t) {
   switch (_type) {
   case (Brownian):
      return _params(1) * (x.array() > -10000).cast <double>().matrix();
   }
}

VectorXd Cestimator::SDE::Model::exact_density(const VectorXd& x0, const VectorXd& xt, double t0, double dt) {
   double mu    = _params(0);
   double sigma = _params(1);
   switch (_type) {
   case (Brownian):
      VectorXd mean = x0.array() + mu * dt;
      VectorXd density(xt.size());

      double sqrt_dt = sqrt(dt);
      for (int i = 0; i < xt.size(); ++i) {
         double diff = xt(i) - mean(i);
         density(i) = (1.0 / (sigma * sqrt_dt * std::sqrt(2 * M_PI))) * std::exp(-0.5 * std::pow(diff / (sigma * sqrt_dt), 2));
      }
      return density;
   }
}

VectorXd Cestimator::SDE::Model::exact_stepping(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
   return x.array() + _params(0) * dt + (dZ * _params(1) * sqrt(dt)).array();
}

VectorXd Cestimator::SDE::Model::dmu_dt(const VectorXd& x, double t) {
   switch (_type) {
   case (Brownian):
      return VectorXd::Zero(x.size());

   default:
      return (1 / eps) * (mu(x, t + eps) - mu(x, t));      //forward diference for time
   }
}

VectorXd Cestimator::SDE::Model::dsigma_dX(const VectorXd& x, double t) {
   switch (_type) {
   case (Brownian):
      return VectorXd::Zero(x.size());

   default:
      return (1 / (2 * eps)) * (sigma(x.array() + eps, t) - sigma(x.array() - eps, t));
   }
}

VectorXd Cestimator::SDE::Model::d2sigma_dX2(const VectorXd& x, double t) {
   switch (_type) {
   case (Brownian):
      return VectorXd::Zero(x.size());

   default:
      return (1 / pow(eps, 2)) * (sigma(x.array() + eps, t) - 2 * sigma(x, t) + sigma(x.array() - eps, t));
   }
}

//and if a model has no clever closed form derivative work, we default to the finite differences
VectorXd Cestimator::SDE::Model::dmu_dX(const VectorXd& x, double t) {
   switch (_type) {
   default:
      return (1 / (2 * eps)) * (mu(x.array() + eps, t) - mu(x.array() - eps, t));         //central differences
   }
}

VectorXd Cestimator::SDE::Model::d2mu_dX2(const VectorXd& x, double t) {
   switch (_type) {
   default:
      return (1 / pow(eps, 2)) * (mu(x.array() + eps, t) - 2 * mu(x, t) + mu(x.array() - eps, t));
   }
}
