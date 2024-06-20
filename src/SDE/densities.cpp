#include "sde.hpp"

inline VectorXd Cestimator::SDE::Density::Exact::operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt) {
   return this->model()->exact_density(x0, xt, t0, dt);
};

inline VectorXd Cestimator::SDE::Density::Euler::operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt) {
   ArrayXd scatter2dt = (this->model()->sigma(x0, t0).array().square()) * 2 * dt;
   ArrayXd location_dt    = x0.array() + this->model()->mu(x0, t0).array() * dt;
   return ((-(xt.array() - location_dt).square() / scatter2dt).exp() / scatter2dt.sqrt() / std::sqrt(M_PI)).matrix();
};

inline VectorXd Cestimator::SDE::Density::Ozaki::operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt) {
   ArrayXd scatter    = this->model()->sigma(x0, t0).array();
   ArrayXd location     = this->model()->mu(x0, t0).array();
   ArrayXd dlocation_dX = this->model()->dmu_dX(x0, t0).array();
   ArrayXd temp   = location * ((dlocation_dX * dt).exp() - 1) / dlocation_dX;

   ArrayXd Mt = x0.array() + temp;
   ArrayXd Kt = (2 / dt) * ((1 + temp / x0.array()).log());
   ArrayXd Vt = scatter * (((Kt * dt).exp() - 1) / Kt).sqrt();

   return ((-0.5 * ((xt.array() - Mt) / Vt).square()).exp() / (std::sqrt(2 * M_PI) * Vt)).matrix();
};

inline VectorXd Cestimator::SDE::Density::ShojiOzaki::operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt) {
   ArrayXd scatter  = this->model()->sigma(x0, t0).array();
   ArrayXd location = this->model()->mu(x0, t0).array();

   ArrayXd Mt = 0.5 * scatter.pow(2) * this->model()->d2mu_dX2(x0, t0).array() + this->model()->dmu_dt(x0, t0).array();
   ArrayXd Lt = this->model()->dmu_dX(x0, t0).array();
   ArrayXd B, A, tmp;
   if ((Lt.cwiseAbs().array() < 1e-5).any()) {
      B = scatter * sqrt(dt);
      A = x0.array() + location * dt + Mt * dt * dt / 2;
   } else {
      B = scatter * (((2 * Lt * dt).exp() - 1)/(2*Lt)).sqrt();
      tmp = (Lt * dt).exp()-1;
      A = x0.array() + location / Lt * tmp + Mt / (Lt * Lt) * (tmp - Lt*dt);
   }
   return ((-0.5 * ((xt.array() - A) / B).pow(2)).exp() / (sqrt(2 * M_PI) * B)).matrix();
   ;
};

inline VectorXd Cestimator::SDE::Density::Elerian::operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt) {
   ArrayXd dscatter_dX = this->model()->dsigma_dX(x0, t0).array();
   if ((dscatter_dX.abs() < 1e-5).any()) {
      return Cestimator::SDE::Density::Euler::operator()(x0, xt, t0, dt);
   }
   ArrayXd scatter  = this->model()->sigma(x0, t0).array();
   ArrayXd location = this->model()->mu(x0, t0).array();

   ArrayXd A = scatter * dscatter_dX * dt * 0.5;
   ArrayXd B = -0.5 * scatter / dscatter_dX + x0.array() + location*dt - A;
   ArrayXd z = (xt.array() - B) / A;
   ArrayXd C = 1 / (dscatter_dX.pow(2) * dt);

   ArrayXd scz = (C*z).sqrt();
   ArrayXd cpz = -0.5 * (C + z);
   ArrayXd ch  = (scz + cpz).exp() + (-scz + cpz).exp();
   return (z.pow(-0.5) * ch / (2 * A.abs() * sqrt(2 * M_PI))).matrix();
};

inline VectorXd Cestimator::SDE::Density::Kessler::operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt) {
   ArrayXd scatter  = this->model()->sigma(x0, t0).array();
   ArrayXd scatter2 = scatter.pow(2);
   ArrayXd location = this->model()->mu(x0, t0).array();

   ArrayXd dscatter_dX   = this->model()->dsigma_dX(x0, t0).array();
   ArrayXd d2scatter_dX2 = this->model()->d2sigma_dX2(x0, t0).array();
   ArrayXd dlocation_dX  = this->model()->dmu_dX(x0, t0).array();

   double   d   = dt * dt / 2;
   ArrayXd E   = x0.array() + location * dt + d * (location * dlocation_dX + 0.5 * scatter2 * d2scatter_dX2);
   ArrayXd tmp = 2 * scatter * dscatter_dX;
   ArrayXd V   = x0.array().pow(2) + (2 * location * x0.array() + scatter2) * dt;
   V += (2 * location * (dlocation_dX * x0.array() + location + scatter * dscatter_dX) +
         scatter2 * (d2scatter_dX2 * x0.array() + 2 * dscatter_dX + tmp + scatter * d2scatter_dX2)) * d - E.pow(2);
   V   = V.abs().sqrt();
   return ((-0.5 * ((xt.array() - E) / V).pow(2)).exp() / (sqrt(2 * M_PI) * V)).matrix();
