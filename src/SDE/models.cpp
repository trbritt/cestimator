#include "sde.hpp"
namespace Cestimator {
namespace SDE {
namespace Model {
class Brownian : public BaseModel {
   Brownian() {
      _is_analytic       = true;
      _simulation_method = "exact";
   };
   VectorXd mu(const VectorXd& x, double t) {
      return _params(0) * (x.array() > -10000).cast <double>().matrix();
   }

   VectorXd sigma(const VectorXd& x, double t) {
      return _params(1) * (x.array() > -10000).cast <double>().matrix();
   }

   VectorXd exact_density(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      double mu    = _params(0);
      double sigma = _params(1);

      VectorXd mean = x0.array() + mu * dt;
      VectorXd density(xt.size());

      double sqrt_dt = sqrt(dt);
      for (int i = 0; i < xt.size(); ++i) {
         double diff = xt(i) - mean(i);
         density(i) = (1.0 / (sigma * sqrt_dt * std::sqrt(2 * M_PI))) * std::exp(-0.5 * std::pow(diff / (sigma * sqrt_dt), 2));
      }
      return density;
   }

   VectorXd exact_stepping(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
      return x.array() + _params(0) * dt + (dZ * _params(1) * sqrt(dt)).array();
   }

   VectorXd dmu_dt(const VectorXd& x, double t) {
      return VectorXd::Zero(x.size());
   }

   VectorXd dsigma_dX(const VectorXd& x, double t) {
      return VectorXd::Zero(x.size());
   }

   VectorXd d2sigma_dX2(const VectorXd& x, double t) {
      return VectorXd::Zero(x.size());
   }
};
}
}
}
