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

   Cestimator::SDE::Model *gbm = new Cestimator::SDE::Model(
      Cestimator::SDE::ModelTypes::Brownian, Cestimator::SDE::PropagatorTypes::Milstein
      );

   gbm->set_params(params);
   std::cout << "hre" << std::endl;

   // now, manage the model by a shared pointer
   std::shared_ptr <Cestimator::SDE::Model> _pgbm(gbm);

   //parameters of simulation
   const double T    = 5;     //num years of sample
   const double freq = 250.0; //num obvs per year
   const double dt   = 1.0 / freq;
   std::cout << dt << std::endl;
   const int M = static_cast <int>(T * freq);

   //create the simualtor
   Cestimator::SDE::Simulator sim(_pgbm, S0, dt, M, 10, 5);

   //run
   MatrixXd paths = sim.simulate_paths();
   std::cout << paths.rows() << " " << paths.cols() << std::endl;
   
   #ifdef __VISUALIZER
   matplot::hold(matplot::on);
   for (int n = 0; n < sim.get_npaths(); ++n) {
      auto l = matplot::plot(Cestimator::Utils::vec2stdvec(VectorXd::LinSpaced(paths.cols(), 0, paths.cols())), Cestimator::Utils::vec2stdvec(paths.row(n)), "-o");
      l->marker_size(5);
      l->marker_face_color("#0000FFF0");
   }
   matplot::show();
   #endif
   Cestimator::Utils::goodbye();
}
