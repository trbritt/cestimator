/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: markowitz.cpp
    license: gpl_v3

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
    MatrixXd data = load_csv<MatrixXd>(fname);
    MatrixXd Y = data(all, all);
    Y.transposeInPlace();

    MatrixXd X = (Y(all, seq(1,last)).array().log()).matrix() - (Y(all, seq(0, last-1)).array().log()).matrix(); //weekly compounded returns

    MatrixXd expected_yr_compounded_returns = X.rowwise().mean() * 52;
    MatrixXd covariance_yr_compounded_returns = Cestimator::Utils::covariance(X) * 52;

    VectorXd current_prices = Y(all, last);

    std::cout << current_prices << " " << std::endl;
    Cestimator::Utils::goodbye();
    return 0;

}