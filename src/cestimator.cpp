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

    MatrixXd data = load_csv<MatrixXd>(fname);

    // std::cout << data(all, last) << std::endl;
    MatrixXd Y = data(all, {dim1,dim2});
    Y.transposeInPlace();

    MatrixXd X = Y(all, seq(1,last)) - Y(all, seq(0, last-1));

    auto res = Estimator::non_parametric(X);
    Vector2d mu_nonpar = std::get<0>(res);
    MatrixXd sigma_nonpar = std::get<1>(res);

    auto res2 = Estimator::shrinkage(X);
    Vector2d mu_shr = std::get<0>(res2);
    MatrixXd sigma_shr = std::get<1>(res2);

    // auto res3 = Estimator::robust(X);

    std::cout << mu_nonpar << " "  << mu_shr << std::endl;
    std::cout << sigma_nonpar << " "  << sigma_shr << std::endl;

    return 0;
}