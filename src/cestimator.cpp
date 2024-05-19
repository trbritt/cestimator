/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: cestimator.cpp
    license: gpl_v3

    this is the implementation for all the things that makes the main
    routine work

*/
#include "cestimator.hpp"

int main(int argc, char *argv[]) {

    if ( argc!=2 ) {
        std::cerr << "Incorrect number of command line arguments passed to programme. Ending" << std::endl;
        return 1;
    } 
    std::string fname = std::string(argv[1]);

    Utils::banner(fname);
    std::vector<Utils::Timer> timers(5);
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

    timers[1].start();
    pprint(Estimator::non_parametric(X), timers[1].name);
    timers[1].stop();

    timers[2].start();
    pprint(Estimator::shrinkage(X), timers[2].name);

    timers[2].stop();

    timers[3].start();
    pprint(Estimator::maximum_likelihood(X), timers[3].name);
    timers[3].stop();

    timers[4].start();
    pprint(Estimator::robust(X), timers[4].name);
    timers[4].stop();

    Utils::goodbye(timers);

    return 0;
}