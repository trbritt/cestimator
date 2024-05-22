#include "cestimator.hpp"

int main(int argc, char *argv[]) {
    if ( argc!=2 ) {
        std::cerr << "Incorrect number of command line arguments passed to programme. Ending" << std::endl;
        return 1;
    } 
    std::string fname = std::string(argv[1]);
    Cestimator::Utils::banner(fname);
    std::vector<Cestimator::Utils::Timer> timers(2);

    timers[0].start();
    MatrixXd data = load_csv<MatrixXd>(fname);

    MatrixXd Y = data(all, all);
    Y.transposeInPlace();
    timers[0].stop();

    Cestimator::Estimator::set_data(Y);
    timers[1].start();
    Cestimator::GMM gmm;
    gmm.name = "GMM";
    gmm.set_data(3);
    gmm.run();
    gmm.print();
    timers[1].stop();

    Cestimator::Utils::goodbye(timers);
}