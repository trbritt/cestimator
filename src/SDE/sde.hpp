#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <optional>
#include <memory>

#define eps    1.0e-6


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
   VectorXd dmu_dX(const VectorXd& x, double t) const {
      return (1 / (2 * eps)) * (mu(x.array() + eps, t) - mu(x.array() - eps, t)); //central differences
   }

   VectorXd dmu_dt(const VectorXd& x, double t) const {
      return (1 / eps) * (mu(x, t + eps) - mu(x, t)); //forward diference for time
   }

   VectorXd d2mu_dX2(const VectorXd& x, double t) const {
      return (1 / pow(eps, 2)) * (mu(x.array() + eps, t) - 2 * mu(x, t) + mu(x.array() - eps, t));
   }

   VectorXd dsigma_dX(const VectorXd& x, double t) const {
      return (1 / (2 * eps)) * (sigma(x.array() + eps, t) - sigma(x.array() - eps, t));
   }

   VectorXd d2sigma_dX2(const VectorXd& x, double t) const {
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

   bool _is_analytic;
   std::string _simulation_method;
   VectorXd _params;
   bool _positive;
};

class BaseDensity {
public:
   BaseDensity(std::shared_ptr <Cestimator::SDE::BaseModel>& model) : _pmodel(model) {
   };
   virtual ~BaseDensity() {
   };
   //overload the call operator
   virtual VectorXd operator()(const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt) = 0;

   template <typename Derived>
   Derived *model() {
      return static_cast <Derived *>(_pmodel.get());
   }

private:
   std::shared_ptr <Cestimator::SDE::BaseModel> _pmodel;  //we need this template to work for all derived models.
   //so we'll create a pointer of the base class and then
   //dynamic cast as necessary later
};
class BasePropagator {
public:
   BasePropagator(std::shared_ptr <Cestimator::SDE::BaseModel>& model) : _pmodel(model) {
   };
   virtual ~BasePropagator() {
   };

   virtual VectorXd next(double t, double dt, const VectorXd& x, const VectorXd& dZ) const = 0;

   VectorXd operator()(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
      return next(t, dt, x, dZ);
   };
   template <typename Derived>
   Derived *model() {
      return static_cast <Derived *>(_pmodel.get());
   }

   std::string id() {
      return _id;
   }

protected:
   std::shared_ptr <Cestimator::SDE::BaseModel> _pmodel;
   std::string _id;
};
class Simulator {
public:
   Simulator(const double S0, const int M, const double dt, const int num_paths, std::shared_ptr <Cestimator::SDE::BasePropagator>& propagator, const std::optional <int> substep = std::nullopt) : _S0(S0), _M(M), _dt(dt), _n_paths(num_paths), _ppropagator(propagator), _substep(substep.value_or(5)) {
   };
   MatrixXd simulate_paths();

private:
   double _S0, _dt;
   int _M, _substep, _n_times, _n_paths;
   MatrixXd _initialise_path(const int num_times);
   MatrixXd _simulate_substep();

protected:
   std::shared_ptr <Cestimator::SDE::BasePropagator> _ppropagator;
};
}
}
