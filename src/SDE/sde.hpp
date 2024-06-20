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
enum ModelTypes { Brownian };
class Model {
public:
   Model(enum ModelTypes type, std::string simulation_method = "milstein")
      : _type(type), _simulation_method(simulation_method), _positive(false) {
      switch (_type) {
      case (Brownian):
         _is_analytic = true;
      }
   };

   //evaluate the drift of the SDE at a point or along a range
   //likewise for diffusion, and since it is all 1D right now,
   //having these as vector types is ok
   VectorXd mu(const VectorXd& x, double t);
   VectorXd sigma(const VectorXd& x, double t);

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
   VectorXd exact_density(const VectorXd& x0, const VectorXd& xt, double t0, double dt);

   VectorXd exact_stepping(double t, double dt, const VectorXd& x, const VectorXd& dZ);

   //and if a model has no clever closed form derivative work, we default to the finite differences
   VectorXd dmu_dX(const VectorXd& x, double t);

   VectorXd dmu_dt(const VectorXd& x, double t);

   VectorXd d2mu_dX2(const VectorXd& x, double t);

   VectorXd dsigma_dX(const VectorXd& x, double t);

   VectorXd d2sigma_dX2(const VectorXd& x, double t);

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
   enum ModelTypes _type;
};

class BaseDensity {
public:
   BaseDensity(Cestimator::SDE::Model& model) : _pmodel(&model) {
   };
   virtual ~BaseDensity() {
   };
   //overload the call operator
   virtual VectorXd operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt) = 0;

   Model *model() {
      return _pmodel;
   }

private:
   Cestimator::SDE::Model *_pmodel;
};
class BasePropagator {
public:
   BasePropagator(Cestimator::SDE::Model& model) : _pmodel(&model) {
   };
   virtual ~BasePropagator() {
   };

   virtual VectorXd next(double t, double dt, const VectorXd& x, const VectorXd& dZ) const = 0;

   VectorXd operator()(double t, double dt, const VectorXd& x, const VectorXd& dZ) {
      return next(t, dt, x, dZ);
   };
   Model *model() {
      return _pmodel;
   }

   std::string id() {
      return _id;
   }

protected:
   Cestimator::SDE::Model *_pmodel;
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

namespace Propagator {
class Exact : public BasePropagator {
   Exact(Cestimator::SDE::Model& model) : BasePropagator(model) {
      _id = "exact";
   }

   VectorXd next(double t, double dt, const VectorXd& x, const VectorXd& dZ);
};

class Euler : public BasePropagator {
   Euler(Cestimator::SDE::Model& model) : BasePropagator(model) {
      _id = "euler";
   }

   VectorXd next(double t, double dt, const VectorXd& x, const VectorXd& dZ);
};

class Milstein : public BasePropagator {
   Milstein(Cestimator::SDE::Model& model) : BasePropagator(model) {
      _id = "milstein";
   }

   VectorXd next(double t, double dt, const VectorXd& x, const VectorXd& dZ);
};
}
namespace Density {
class Exact : public BaseDensity {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt);
};

class AitSahalia : public BaseDensity {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt);
};

class Euler : public BaseDensity {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt);
};

class Ozaki : public BaseDensity {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt);
};
class ShojiOzaki : public BaseDensity {
public:
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt);
};

class Elerian : public BaseDensity {
public:
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt);
};

class Kessler : public BaseDensity {
   VectorXd operator()(const VectorXd& x0, const VectorXd& xt, double t0, double dt);
};
}
}
}
