#include "sde.hpp"

VectorXd Cestimator::SDE::Model::next(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
   VectorXd xp = VectorXd::Zero(x.size());
   switch (_propagatortype) {
   case (Exact):
   {
      xp = exact_stepping(t, dt, x, dZ);
      break;
   }

   case (Euler):
   {
      xp = x.array() + mu(x, t).array() * dt + sigma(x, t).array() * sqrt(dt) * dZ.array();
      xp = positivity() ? xp.array().max(0.0) : xp;
      break;
   }

   case (Milstein):
   {
      xp = x.array() + mu(x, t).array() * dt +
           sigma(x, t).array() * sqrt(dt) * dZ.array() +
           0.5 * sigma(x, t).array() * dsigma_dX(x, t).array() * (dZ.cwiseProduct(dZ).array() - 1) * dt;
      xp = positivity() ? xp.array().max(0.0) : xp;
      break;
   }
   }
   return xp;
}
