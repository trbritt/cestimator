#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <optional>

double eps = 1.0e-6;

using namespace Eigen;

/*
 * Create abstract Model class onto which we will inherit the specific models
 */
namespace Cestimator {
namespace SDE {
class BaseModel {
public:
   BaseModel(const std::optional <bool> is_analytic = std::nullopt, std::string simulation_method = "milstein")
      : _is_analytic(is_analytic), _simulation_method(simulation_method), _positive(false) {
   };
   virtual ~BaseModel() = default; //virtual destructor to ensure deconstruct derived children correctly

   //evaluate the drift of the SDE at a point or along a range
   //likewise for diffusion, and since it is all 1D right now,
   //having these as vector types is ok
   virtual VectorXd mu(const VectorXd& x, double t) const    = 0;
   virtual VectorXd sigma(const VectorXd& x, double t) const = 0;

   VectorXd get_params() const {
      return _params;
   };
   bool has_analytic_density() const {
      return _is_analytic;
   }

   // these are all the density estimators we will override for each model
   //because we can either evaluate these at points or vectors,
   // we'll only need to do this for the exact density of analytic models
   //and therefore the Hermitian expansion of those forms
   virtual VectorXd exact_density(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) = 0;

   virtual VectorXd ait_sahalia_density(const VectorXd&  x0, const VectorXd&  xt, const VectorXd&  t0, double dt) = 0;

   virtual VectorXd exact_stepping(double t, double dt, const VectorXd& x, const VectorXd& dZ) = 0;

   //and if a model has no clever closed form derivative work, we default to the finite differences
   VectorXd dmu_dX(VectorXd& x, double t) const {
      return (1 / (2 * eps)) * (mu(x.array() + eps, t) - mu(x.array() - eps, t)); //central differences
   }

   VectorXd dmu_dt(VectorXd& x, double t) const {
      return (1 / eps) * (mu(x, t + eps) - mu(x, t)); //forward diference for time
   }

   VectorXd d2mu_dX2(VectorXd& x, double t) const {
      return (1 / pow(eps, 2)) * (mu(x.array() + eps, t) - 2 * mu(x, t) + mu(x.array() - eps, t));
   }

   VectorXd dsigma_dX(VectorXd& x, double t) const {
      return (1 / (2 * eps)) * (sigma(x.array() + eps, t) - sigma(x.array() - eps, t));
   }

   VectorXd d2sigma_dX2(VectorXd& x, double t) const {
      return (1 / pow(eps, 2)) * (sigma(x.array() + eps, t) - 2 * sigma(x, t) + sigma(x.array() - eps, t));
   }

   void set_params(const VectorXd& p) {
      _positive = set_positivity(p);
      _params   = p;
   };

   bool positivity() const {
      return _positive;
   }

   std::string simulation_method() const {
      return _simulation_method;
   }

protected:
   virtual bool set_positivity(const VectorXd& p) const {
      return false;
   };

private:
   bool _is_analytic;
   std::string _simulation_method;
   VectorXd _params;
   bool _positive;
};

class Density {
public:
   Density(Cestimator::SDE::BaseModel& model) : _pmodel(&model) {
   };
   virtual ~Density() {
   };
   //overload the call operator
   virtual VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) = 0;

   template <typename Derived>
   Derived *model() {
      return dynamic_cast <const Derived *>(*_pmodel);
   }

private:
   Cestimator::SDE::BaseModel *_pmodel;  //we need this template to work for all derived models.
                                         //so we'll create a pointer of the base class and then
                                         //dynamic cast as necessary later
};
}
}
