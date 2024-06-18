#include "sde.hpp"

namespace Cestimator {
namespace SDE {
template <typename Derived>
class Exact : public Density {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      return this->model <Derived>()->exact_density(x0, xt, t0, dt);
   }
};

template <typename Derived>
class AitSahalia : public Density {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      return this->model <Derived>()->ait_sahalia_density(x0, xt, t0, dt);
   }
};

template <typename Derived>
class Euler : public Density {
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
class Ozaki : public Density {
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
template <typename Derived>
class ShojiOzaki : public Density {
public:
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      VectorXd scatter  = this->model <Derived>()->sigma(x0, t0);
      VectorXd location = this->model <Derived>()->mu(x0, t0);

      VectorXd Mt = 0.5 * scatter.cwiseProduct(scatter) * this->model <Derived>()->d2mu_dX2(x0, t0) + this->model <Derived>()->dmu_dt(x0, t0);
      VectorXd Lt = this->model <Derived>()->dmu_dX(x0, t0);
      VectorXd B, A, tmp;
      if ((Lt.cwiseAbs().array() < 1e-5).any()) {
         B = scatter * sqrt(dt);
         A = x0 + location * dt + Mt * dt * dt / 2;
      } else {
         B = (2 * Lt * dt).unaryExpr([](double x) {
                  return exp(x) - 1;
               }).cwiseQuotient(2 * Lt);
         tmp = (Lt * dt).unaryExpr([](double x) {
                  return exp(x);
               }) - 1;
         A = x0 + location.cwiseQuotient(Lt).cwiseProduct(tmp) + Mt.cwiseQuotient(Lt.cwiseProduct(Lt)).cwiseProduct(tmp - Lt * dt);
      }
      tmp = (xt - A).cwiseQuotient(B);
      VectorXd numerator = tmp.unaryExpr([](double x) {
               return exp(-0.5 * x * x);
            });
      VectorXd denominator = B * sqrt(2 * M_PI);
      return numerator.cwiseQuotient(denominator);
   }
};

template <typename Derived>
class Elerian : public Density {
public:
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      VectorXd dscatter_dX = this->model <Derived>()->dsigma_dX(x0, t0);
      if ((dscatter_dX.cwiseAbs().array() < 1e-5).any()) {
         return Cestimator::SDE::Euler <Derived>::operator()(x0, xt, t0, dt);
      }
      VectorXd scatter  = this->model <Derived>()->sigma(x0, t0);
      VectorXd location = this->model <Derived>()->mu(x0, t0);

      VectorXd A = scatter.cwiseProduct(dscatter_dX) * dt / 2;
      VectorXd B = -0.5 * scatter.cwiseQuotient(dscatter_dX) + x0 + location * dt - A;
      VectorXd z = (xt - B).cwiseQuotient(A);
      VectorXd C = (dscatter_dX.cwiseProduct(dscatter_dX) * dt).cwiseInverse();

      VectorXd scz = C.cwiseProduct(z).cwiseSqrt();
      VectorXd cpz = -0.5 * (C + z);
      VectorXd ch  = (scz + cpz).unaryExpr([](double x) {
               return exp(x);
            }) + (-scz + cpz).unaryExpr([](double x) {
               return exp(x);
            });
      return z.cwiseSqrt().cwiseInverse().cwiseProduct(ch).cwiseQuotient(2 * A.cwiseAbs() * sqrt(2 * M_PI));
   }
};

template <typename Derived>
class Kessler : public Density {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) {
      VectorXd scatter  = this->model <Derived>()->sigma(x0, t0);
      VectorXd scatter2 = scatter.cwiseProduct(scatter);
      VectorXd location = this->model <Derived>()->mu(x0, t0);

      VectorXd dscatter_dX   = this->model <Derived>()->dsigma_dX(x0, t0);
      VectorXd d2scatter_dX2 = this->model <Derived>()->d2sigma_dX2(x0, t0);
      VectorXd dlocation_dX  = this->model <Derived>()->dmu_dX(x0, t0);

      double   d   = dt * dt / 2;
      VectorXd E   = x0 + location * dt + d * (location.cwiseProduct(dlocation_dX) + 0.5 * scatter2.cwiseProduct(d2scatter_dX2));
      VectorXd tmp = 2 * scatter.cwiseProduct(dscatter_dX);
      VectorXd V   = x0.cwiseProduct(x0) + (2 * location.cwiseProduct(x0) + scatter2) * dt;
      V += (2 * location * (dlocation_dX.cwiseProduct(x0) + location + scatter.cwiseProduct(dscatter_dX)) +
            scatter2.cwiseProduct(d2scatter_dX2.cwiseProduct(x0) + 2 * dscatter_dX + tmp + scatter.cwiseProduct(d2scatter_dX2))) * d - E.cwiseProduct(E);
      V   = V.cwiseAbs().cwiseSqrt();
      tmp = (xt - E).cwiseQuotient(V);
      return tmp.unaryExpr([](double x) {
               return exp(-0.5 * x * x);
            }).cwiseQuotient(V * sqrt(2 * M_PI));
   }
};
}
}
