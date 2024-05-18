#include "cestimator.hpp"
#include "support.hpp"

int main(int argc, char *argv[]) {

    if ( argc!=4 ) {
        std::cerr << "Incorrect number of command line arguments passed to programme. Ending" << std::endl;
        return 1;
    } 
    std::string fname = std::string(argv[1]);
    int dim1 = std::stoi(argv[2]);
    int dim2 = std::stoi(argv[3]);

    Eigen::MatrixXd data = load_csv<Eigen::MatrixXd>(fname);

    // std::cout << data(all, last) << std::endl;
    Eigen::MatrixXd Y = data(all, {dim1,dim2});
    Y.transposeInPlace();

    Eigen::MatrixXd X = Y(all, seq(1,last)) - Y(all, seq(0, last-1));

    auto res = non_parametric_estimation(X);
    Vector2d mu_nonpar = std::get<0>(res);
    MatrixXd sigma_nonpar = std::get<1>(res);

    auto res2 = shrinkage_estimation(X);
    Vector2d mu_shr = std::get<0>(res2);
    MatrixXd sigma_shr = std::get<1>(res2);

    // auto res3 = hubert_M(X);

    std::cout << mu_nonpar << " "  << mu_shr << std::endl;
    std::cout << sigma_nonpar << " "  << sigma_shr << std::endl;

    return 0;
}