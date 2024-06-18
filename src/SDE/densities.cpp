#include "sde.hpp"

namespace Cestimator{
    namespace SDE{
        template<typename Derived>
        class ExactDensity: public Density {
            VectorXd operator() (const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt){
                return this->model<Derived>()->exact_density(x0, xt, t0, dt);
            }
        };

        template<typename Derived>
        class EulerDensity: public Density {
            VectorXd operator() (const VectorXd& x0, const VectorXd& xt, const VectorXd& t0, double dt){
                VectorXd diffusion_contribution = this->model<Derived>()->sigma(x0, t0).cwiseProduct(this->model<Derived>()->sigma(x0, t0)) * 2 * dt;
                VectorXd drift_contribution = x0 + this->model<Derived>()->mu(x0, t0)*dt;
                
                VectorXd numerator = (xt - drift_contribution).cwiseProduct(xt-drift_contribution).cwiseQuotient(diffusion_contribution).unaryExpr([](double x) {
                    return exp(-x);
                });
                VectorXd denominator = diffusion_contribution.cwiseSqrt() * sqrt(M_PI);
                return numerator.cwiseQuotient(denominator);
            }
        };

    }
}