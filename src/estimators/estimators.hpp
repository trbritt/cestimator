/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: estimators.hpp
    license: gpl_v3

    this is the header for all of the estimators currently implemented in the programme

*/
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>
#include <iterator>
#include "../utils/utils.hpp"

using namespace Eigen;

namespace Cestimator{

    class Estimator {
            
        protected:
            static MatrixXd data,sigma_no_par;
            static VectorXd mu_no_par;
            static int N, T;

        public:

            VectorXd mu;
            MatrixXd sigma;

            std::string name;
            Estimator(){};
            static inline void set_data(MatrixXd& arr){

                data = arr;

                N = data.rows();
                T = data.cols();

                mu_no_par = Utils::mean(data) ;
                sigma_no_par = Utils::covariance(data);
            }
            inline void print(){
                size_t num_spaces = std::max(0, static_cast<int>(15-name.length()));
                std::string padding(num_spaces, ' ');

                int ndims = mu.size();
                std::cout << Cestimator::Utils::colors::OKCYAN << name + padding << "\t\t" << "μ" << "\t\t" << "Σ" << Cestimator::Utils::colors::ENDC << std::endl;

                for (int i=0; i<ndims; ++i){
                    std::cout << "\t\t" <<  mu(i) << "\t";
                    for (int j=0; j<ndims; ++j){
                        std::cout << sigma(i,j) << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "\n" << std::endl;
            };
            virtual int run() noexcept = 0 ; //pure virtual initialization, no one can inheret unless they impl run. label as noexcept since this should not throw any issues
    };

    typedef std::tuple<VectorXd, MatrixXd> Result;

    class non_parametric: public Estimator{
        public:
            int run() noexcept override; 
    };

    class shrinkage: public Estimator{
        public:
            int run() noexcept override;
    };

    class maximum_likelihood: public Estimator{
        public:
            int run() noexcept override;
    };

    class robust: public Estimator{
        public:
            int run() noexcept override;
    };

    class GMM: public Estimator {
        public:
            int n_features, n_term=-1;
            MatrixXd mu, regularization;
            std::vector<MatrixXd> sigma;
            inline void set_data(int nf=3){
                n_features = nf;
            }
            int run() noexcept override;
            VectorXd predict(VectorXd& arr);
            void print();
    };
        
}