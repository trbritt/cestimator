#include "sde.hpp"

/*
 * Brownian:
 *  params: mu, sigma
 *  SDE: dX_t = mu(X,t)dt + sigma(X,t)dW_t
 *                 =             =
 *                 mu           sigma > 0
 * CEV:
 *  params: kappa, mu, sigma, gamma
 *  SDE: dX_t = mu(X,t)dt + sigma(X,t)dW_t
 *                  =           =
 *              kappa(mu-X)  sigma*X^gamma
 *
 *
 */
VectorXd Cestimator::SDE::Model::mu(const VectorXd& x, double t) {
   VectorXd retval = VectorXd::Zero(x.size());
   switch (_modeltype) {
   case (Brownian):
   {
      retval = _params(0) * (x.array() > -10000).cast <double>().matrix();
      break;
   }

   case (CEV):
   {
      retval = _params(0) * (_params(1) - x.array()).matrix();
      break;
   }
   }
   return retval;
}

VectorXd Cestimator::SDE::Model::sigma(const VectorXd& x, double t) {
   VectorXd retval = VectorXd::Zero(x.size());

   switch (_modeltype) {
   case (Brownian):
   {
      retval = _params(1) * (x.array() > -10000).cast <double>().matrix();
      break;
   }

   case (CEV):
   {
      retval = _params(2) * x.array().pow(_params(3)).matrix();
      break;
   }
   }
   return retval;
}

VectorXd Cestimator::SDE::Model::exact_density(const VectorXd& x0, const VectorXd& xt, double t0, double dt) {
   double   mu      = _params(0);
   double   sigma   = _params(1);
   VectorXd density = VectorXd::Zero(xt.size());
   VectorXd mean;
   switch (_modeltype) {
   case (Brownian):
   {
      mean = x0.array() + mu * dt;
      double sqrt_dt = sqrt(dt);
      for (int i = 0; i < xt.size(); ++i) {
         double diff = xt(i) - mean(i);
         density(i) = (1.0 / (sigma * sqrt_dt * std::sqrt(2 * M_PI))) * std::exp(-0.5 * std::pow(diff / (sigma * sqrt_dt), 2));
      }
      break;
   }

   default:
      std::cerr << "Exact density not available for this model" << std::endl;
      break;
   }
   return density;
}

VectorXd Cestimator::SDE::Model::exact_stepping(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
   VectorXd retval = VectorXd::Zero(x.size());
   switch (_modeltype) {
   case (Brownian):
   {
      retval = x.array() + _params(0) * dt + (dZ * _params(1) * sqrt(dt)).array();
      break;
   }

   default:
      std::cerr << "Exact stepping cannot be used for this type of model" << std::endl;
      break;
   }
   return retval;
}

VectorXd Cestimator::SDE::Model::dmu_dt(const VectorXd& x, double t) {
   VectorXd retval = VectorXd::Zero(x.size());

   switch (_modeltype) {
   case (Brownian):
      break;

   default:
   {
      retval = (1 / eps) * (mu(x, t + eps) - mu(x, t));        //forward diference for time
      break;
   }
   }
   return retval;
}

VectorXd Cestimator::SDE::Model::dsigma_dX(const VectorXd& x, double t) {
   VectorXd retval = VectorXd::Zero(x.size());

   switch (_modeltype) {
   case (Brownian):
      break;

   default:
   {
      retval = (1 / (2 * eps)) * (sigma(x.array() + eps, t) - sigma(x.array() - eps, t));
      break;
   }
   }
   return retval;
}

VectorXd Cestimator::SDE::Model::d2sigma_dX2(const VectorXd& x, double t) {
   VectorXd retval = VectorXd::Zero(x.size());

   switch (_modeltype) {
   case (Brownian):
      break;

   default:
   {
      retval = (1 / pow(eps, 2)) * (sigma(x.array() + eps, t) - 2 * sigma(x, t) + sigma(x.array() - eps, t));
      break;
   }
   }
   return retval;
}

//and if a model has no clever closed form derivative work, we default to the finite differences
VectorXd Cestimator::SDE::Model::dmu_dX(const VectorXd& x, double t) {
   switch (_modeltype) {
   default:
      return (1 / (2 * eps)) * (mu(x.array() + eps, t) - mu(x.array() - eps, t));         //central differences
   }
}

VectorXd Cestimator::SDE::Model::d2mu_dX2(const VectorXd& x, double t) {
   switch (_modeltype) {
   default:
      return (1 / pow(eps, 2)) * (mu(x.array() + eps, t) - 2 * mu(x, t) + mu(x.array() - eps, t));
   }
}
