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

    if ( argc!=4 ) {
        std::cerr << "Incorrect number of command line arguments passed to programme. Ending" << std::endl;
        return 1;
    } 
    std::string fname = std::string(argv[1]);
    int dim1 = std::stoi(argv[2]);
    int dim2 = std::stoi(argv[3]);

    Utils::banner(fname, dim1, dim2);
    std::vector<Utils::Timer> timers(5);
    timers[0].name = "loading";
    timers[1].name = "nonparametric";
    timers[2].name = "shrinkage";
    timers[3].name = "MLE";
    timers[4].name = "robust";

    timers[0].start();
    MatrixXd data = load_csv<MatrixXd>(fname);

    // std::cout << data(all, last) << std::endl;
    MatrixXd Y = data(all, {dim1,dim2});
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
    // auto res3 = Estimator::robust(X);
    timers[3].stop();

    timers[4].start();
    timers[4].stop();

    Utils::goodbye(timers);

    return 0;
}