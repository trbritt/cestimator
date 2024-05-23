/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: basic.cpp
    license: gpl_v3

    this is a basic example illustrating the different estimators

*/
#include "cestimator.hpp"

int main(int argc, char *argv[]) {

    if ( argc!=2 ) {
        std::cerr << "Incorrect number of command line arguments passed to programme. Ending" << std::endl;
        return 1;
    } 
    
    Eigen::initParallel();
    Eigen::setNbThreads(6);
    
    std::string fname = std::string(argv[1]);

    Cestimator::Utils::banner(fname);
    std::vector<Cestimator::Utils::Timer> timers(5);
    timers[0].name = "loading";
    timers[1].name = "nonparametric";
    timers[2].name = "shrinkage";
    timers[3].name = "MLE";
    timers[4].name = "robust";

    timers[0].start();
    MatrixXd data = load_csv<MatrixXd>(fname);

    MatrixXd Y = data(all, all);
    Y.transposeInPlace();

    MatrixXd X = Y(all, seq(1,last)) - Y(all, seq(0, last-1));
    timers[0].stop();

    Cestimator::Estimator::set_data(X);
    Cestimator::non_parametric np;
    Cestimator::shrinkage shr;
    Cestimator::maximum_likelihood mle;
    Cestimator::robust rb;

    std::vector<Cestimator::Estimator*> estimators = {&np, &shr, &mle, &rb};
    int i=0;
    for ( auto est : estimators){
        timers[i+1].start();
        est->name = timers[i+1].name;
        est->run();
        est->print();
        timers[i+1].stop();
        ++i;
    }

    #ifdef __VISUALIZER
    int dim1=0, dim2=2;
    Cestimator::Utils::Visualizer viz;
    viz.scatter2d(X, dim1, dim2);
    viz.ellipse(rb.mu, rb.sigma, dim1, dim2);
    viz.ellipse(mle.mu, mle.sigma, dim1, dim2);
    viz.ellipse(shr.mu, shr.sigma, dim1, dim2);
    viz.ellipse(np.mu, np.sigma, dim1, dim2);
    viz.gr->Run();
    #endif

    Cestimator::Utils::goodbye(timers);
    return 0;
}
