#include "sde.hpp"

namespace Cestimator {
namespace SDE {
template <typename Derived>
class ExactDensity : public Density {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      return this->model <Derived>()->exact_density(x0, xt, t0, dt);
   }
};

template <typename Derived>
class AitSahaliaDensity : public Density {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      return this->model <Derived>()->ait_sahalia_density(x0, xt, t0, dt);
   }
};

template <typename Derived>
class EulerDensity : public Density {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      VectorXd diffusion_contribution = this->model <Derived>()->sigma(x0, t0).cwiseProduct(this->model <Derived>()->sigma(x0, t0)) * 2 * dt;
      VectorXd drift_contribution     = x0 + this->model <Derived>()->mu(x0, t0) * dt;

      VectorXd numerator = (xt - drift_contribution).cwiseProduct(xt - drift_contribution).cwiseQuotient(diffusion_contribution).unaryExpr([](double x) {
               return exp(-x);
            });
      VectorXd denominator = diffusion_contribution.cwiseSqrt() * sqrt(M_PI);
      return numerator.cwiseQuotient(denominator);
   }
};

template <typename Derived>
class OzakiDensity : public Density {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      VectorXd scatter      = this->model <Derived>()->sigma(x0, t0);
      VectorXd location     = this->model <Derived>()->mu(x0, t0);
      VectorXd dlocation_dX = this->model <Derived>()->dmu_dX(x0, t0);

      VectorXd tmp = ((dlocation_dX * dt).unaryExpr([](double x) {
               return exp(x);
            }) - 1) / dlocation_dX;

      VectorXd Mt = x0 + tmp;
      VectorXd Kt = (2 / dt) * (1 + tmp.cwiseQuotient(x0).array()).unaryExpr([](double x) {
               return log(x);
            });
      VectorXd Vt = ((Kt * dt).unaryExpr([](double x) {
               return exp(x);
            }) - 1).cwiseQuotient(Kt).cwiseSqrt() * scatter;

      VectorXd numerator = (xt - Mt).cwiseQuotient(Vt).cwiseProduct((xt - Mt).cwiseQuotient(Vt)).unaryExpr([](double x) {
               return exp(-0.5 * x);
            });
      VectorXd denominator = Vt * sqrt(M_PI * 2);
      return numerator.cwiseQuotient(denominator);
   }
};
}
}
