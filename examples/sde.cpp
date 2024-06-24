#include "cestimator.hpp"

int main(int argc, char *argv[]) {
   Eigen::initParallel();
   Eigen::setNbThreads(6);

   Cestimator::Utils::banner("");

   const double S0    = 0.4;
   const double mu    = 2;
   const double sigma = 0.5;

   VectorXd params(2);
   params << mu, sigma;

   try {
      Cestimator::SDE::Model *gbm = new Cestimator::SDE::Model(
         Cestimator::SDE::ModelTypes::Brownian, "milstein"
         );

      gbm->set_params(params);
      std::cout << "hre" << std::endl;

      // now, manage the model by a shared pointer
      std::shared_ptr <Cestimator::SDE::Model> _pgbm(gbm);

      //create propagator
      Cestimator::SDE::Propagator *propogator = new Cestimator::SDE::Propagator(Cestimator::SDE::PropagatorTypes::Exact, _pgbm);

      //and shared ptr to it
      std::shared_ptr <Cestimator::SDE::Propagator> _peprop(propogator);

      //parameters of simulation
      const double T    = 5;   //num years of sample
      const int    freq = 250; //num obvs per year
      const double dt   = 1 / freq;
      const int    M    = static_cast <int>(T * freq);

      //create the simualtor
      Cestimator::SDE::Simulator sim(S0, M, dt, 10, _peprop, 5);

      //run
      MatrixXd paths = sim.simulate_paths();
      std::cout << paths.rows() << " " << paths.cols() << std::endl;
   }
   catch (const std::bad_weak_ptr& bwp) {
      std::cerr << bwp.what();
      exit(-1);
   }

   Cestimator::Utils::goodbye();
}
