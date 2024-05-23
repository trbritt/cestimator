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

    timers[1].start();
    Cestimator::non_parametric np;
    np.name = timers[1].name;
    np.run();
    np.print();
    timers[1].stop();

    timers[2].start();
    Cestimator::shrinkage shr;
    shr.name = timers[2].name;
    shr.run();
    shr.print();
    timers[2].stop();

    timers[3].start();
    Cestimator::maximum_likelihood mle;
    mle.name = timers[3].name;
    mle.run();
    mle.print();
    timers[3].stop();

    timers[4].start();
    Cestimator::robust rb;
    rb.name = timers[4].name;
    rb.run();
    rb.print();
    timers[4].stop();

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
