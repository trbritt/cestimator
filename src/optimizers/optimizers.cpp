#include "optimizers.hpp"
#include <random>
std::random_device rd;     // Only used once to initialise (seed) engine
std::mt19937 rng(rd());  
std::uniform_int_distribution<int> uni(0,1); // Guaranteed unbiased

namespace Cestimator {

    VectorXd Optimizer::generate_weights(VectorXd returns, MatrixXd covariance){
        int n_assets = returns.size();
        VectorXd weights = VectorXd::Zero(n_assets);
        switch ( static_cast<int>(classifier))
        {
        case 0: //random weights
            for (int i = 0; i<n_assets; ++i){
                weights(i) = uni(rng);
            }
        case 1:
            break;
        }
    }

}