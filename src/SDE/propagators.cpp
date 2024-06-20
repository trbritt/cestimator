#include "sde.hpp"

namespace Cestimator {
namespace SDE {
namespace Propagator {
template <typename Derived>
class Exact : public BasePropagator {
   Exact(std::shared_ptr <Cestimator::SDE::BaseModel>& model) : BasePropagator(model) {
      _id = "exact";
   }

   VectorXd next(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
      return this->model <Derived>()->exact_stepping(t, dt, x, dZ);
   }
};

template <typename Derived>
class Euler : public BasePropagator {
   Euler(std::shared_ptr <Cestimator::SDE::BaseModel>& model) : BasePropagator(model) {
      _id = "euler";
   }

   VectorXd next(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
      VectorXd xp = x.array + this->model <Derived>()->mu(x, t) * dt + this->model <Derived>()->sigma(x, t) * sqrt(dt) * dZ;
      return this->model <Derived>()->_positive ? xp.array().max(0.0) : xp;
   }
};

template <typename Derived>
class Milstein : public BasePropagator {
   Milstein(std::shared_ptr <Cestimator::SDE::BaseModel>& model) : BasePropagator(model) {
      _id = "milstein";
   }

   VectorXd next(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
      VectorXd xp = x.array + this->model <Derived>()->mu(x, t) * dt +
                    this->model <Derived>()->sigma(x, t) * sqrt(dt) * dZ +
                    0.5 * this->model <Derived>()->sigma(x, t) * this->model <Derived>()->dsigma_dX(x, t) * (dZ.cwiseProduct(dZ).array() - 1) * dt;
      return this->model <Derived>()->_positive ? xp.array().max(0.0) : xp;
   }
};
}
}
}
