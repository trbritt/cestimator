/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: cestimator.cpp
    license: gpl_v3

    this is the implementation for all the things that makes the main
    routine work

*/
#include "cestimator.hpp"

MatrixXd X, Y; // put in global name space so dont have to be captured by lambdas

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

    Y = data(all, all);
    Y.transposeInPlace();

    X = Y(all, seq(1,last)) - Y(all, seq(0, last-1));
    timers[0].stop();

    timers[1].start();
    pprint(Cestimator::non_parametric(X), timers[1].name);
    timers[1].stop();

    timers[2].start();
    pprint(Cestimator::shrinkage(X), timers[2].name);

    timers[2].stop();

    timers[3].start();
    pprint(Cestimator::maximum_likelihood(X), timers[3].name);
    timers[3].stop();

    timers[4].start();
    pprint(Cestimator::robust(X), timers[4].name);
    timers[4].stop();

    
    mglFLTK gr([](auto a){return Cestimator::Utils::Visualizers::scatter3d(a, X);}, "DataScatter");
    int success = gr.Run();
    Cestimator::Utils::goodbye(timers);
    return 0;
}