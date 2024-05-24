#include "../utils/utils.hpp"
#include <optional>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

namespace Cestimator{

    class Evaluator {
        
        public:
            enum Distribution {Normal, Uniform, StudentT};
            static MatrixXd arr, sigma;
            static VectorXd mu;
            static Distribution dist;
            const static int N,T;

            Evaluator(enum Distribution d){
                dist = d;
            };
            void set_data(MatrixXd x, VectorXd m, MatrixXd s){
                arr = x;
                mu = m;
                sigma = s;
                const int N = arr.rows();
                const int T = arr.cols();
            };
            MatrixXd generator(const std::optional<int> nu = std::nullopt){
                MatrixXd inv_sigma = sigma.llt().solve(MatrixXd::Identity(N, N));
                MatrixXd ma2 = ((arr.colwise()-mu).transpose() * inv_sigma * (arr.colwise()-mu)).rowwise().sum();

                switch (dist) {
                    case Uniform:
                        return (ma2.array() <= 1 && ma2.array() >= 0).select(
                            0.0,
                            MatrixXd::Ones(N,T)*gamma(static_cast<double>(N)/2+1) * pow(M_PI, -static_cast<double>(N)/2)
                        );
                    case Normal:
                        return (-0.5*ma2).exp();
                    case StudentT:
                        double dof = static_cast<double>(nu.value_or(1));
                        return (1+ma2.array()/dof).matrix().pow(-0.5*(dof+static_cast<double>(N))) * gamma(0.5*(dof+static_cast<double>(N))) / (gamma(0.5*dof) * pow(dof*M_PI, 0.5*N));
                    default:
                        return MatrixXd::Ones(N,T);
                };
            };
    };

}