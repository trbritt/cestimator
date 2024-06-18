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

class BaseModel {
public:
   BaseModel(const std::optional <bool> is_analytic = std::nullopt, std::string simulation_method = "milstein")
      : _is_analytic(is_analytic), _simulation_method(simulation_method), _positive(false) {
   };
   virtual ~BaseModel() = default; //virtual destructor to ensure deconstruct derived children correctly

   //evaluate the drift of the SDE at a point or along a range
   //likewise for diffusion, and since it is all 1D right now,
   //having these as vector types is ok
   virtual double mu(double x, double t) const = 0;
   virtual VectorXd drift(const VectorXd& x, double t) const = 0;

   virtual double sigma(double x, double t) const      = 0;
   virtual VectorXd sigma(const VectorXd& x, double t) = 0;

   VectorXd get_params() const {
      return _params;
   };
   bool has_analytic_density() const {
      return _is_analytic;
   }

   // these are all the density estimators we will override for each model
   virtual double exact_density(double x0, double xt, double t0, double dt)           = 0;
   virtual double ait_sahalia_density(double x0, double xt, double t0, double dt)     = 0;
   virtual double exact_stepping(double t, double dt, double x, double dZ)            = 0;
   virtual VectorXd exact_stepping(double t, double dt, const VectorXd& x, double dZ) = 0;

   //and if a model has no clever closed form derivative work, we default to the finite differences
   double dmu_dX(double x, double t) const {
      return (1 / (2 * eps)) * (mu(x + eps, t) - mu(x - eps, t)); //central differences
   }

   double dmu_dt(double x, double t) const {
      return (1 / eps) * (mu(x, t + eps) - mu(x, t)); //forward diference for time
   }

   double d2mu_dX2(double x, double t) const {
      return (1 / pow(eps, 2)) * (mu(x + eps, t) - 2 * mu(x, t) + mu(x - eps, t));
   }

   double dsigma_dX(double x, double t) const {
      return (1 / (2 * eps)) * (sigma(x + eps, t) - sigma(x - eps, t));
   }

   double d2sigma_dX2(double x, double t) const {
      return (1 / pow(eps, 2)) * (sigma(x + eps, t) - 2 * sigma(x, t) + sigma(x - eps, t));
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
