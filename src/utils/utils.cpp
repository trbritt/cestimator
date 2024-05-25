/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: utils.cpp
    license: gpl_v3

    this is the implementation for all utilities of the project,
    namely anything not theoretically driving the main results,
    but necessary to get there

*/

#include "utils.hpp"

using namespace Eigen;

namespace Cestimator{
    
    namespace Utils {

        const std::string colors::HEADER = "\033[95m";
        const std::string colors::OKBLUE = "\033[94m";
        const std::string colors::OKCYAN = "\033[96m";
        const std::string colors::OKGREEN = "\033[92m";
        const std::string colors::WARNING = "\033[93m";
        const std::string colors::FAIL = "\033[91m";
        const std::string colors::ENDC = "\033[0m";
        const std::string colors::BOLD = "\033[1m";
        const std::string colors::UNDERLINE = "\033[4m";

        VectorXd mean(MatrixXd& x){
            return x.rowwise().mean();
        }

        MatrixXd covariance(MatrixXd& x){
            MatrixXd tmp = x.transpose();
            MatrixXd centered = tmp.rowwise() - tmp.colwise().mean();
            return (tmp.adjoint() * centered) / (double(tmp.rows() - 1));
        }

        VectorXd outlier_cutoff(VectorXd& d, double d0){
            const double b_2 = 1.25;
            const double b_2_2 = b_2 * b_2;

            VectorXd omega = (d.array() <= d0).select(
                1.0,
                (-0.5 * ((d.array() - d0).square() / b_2_2)).exp()
            );
            return omega;
        }

        std::vector<double> matrix2vector(MatrixXd& x){
            MatrixXd xx = x;
            int N, M;
            N = x.rows();
            M = x.cols();

            std::vector<double> data;
            xx.resize(1, N*M);
            for(int i=0; i<N*M; ++i){
                data.push_back(xx(i));
            }
            return data;
        }
    }
}