#include "optimizers.hpp"
#include "eiquadprog.hpp"
#include <random>
std::random_device rd;     // Only used once to initialise (seed) engine
std::mt19937 rng(rd());  
std::uniform_int_distribution<int> uni(0,1); // Guaranteed unbiased

namespace Cestimator {

    MatrixXd Optimizer::generate_frontier(VectorXd returns, MatrixXd covariance, int n_portfolios){
        int n_assets = returns.size();
        MatrixXd weights = MatrixXd::Zero(n_portfolios, n_assets);
        switch ( static_cast<int>(classifier))
        {
        case 0: //random weights
            for (int p=0; p<n_portfolios; ++p){
                for (int i = 0; i<n_assets; ++i){
                    weights(p, i) = uni(rng);
                }
                profile.push_back(
                    std::make_pair(
                        weights.row(p).dot(returns),
                        sqrt(weights.row(p)*covariance*weights.row(p).transpose())
                    )
                );
            }
            
        case 1:
            MatrixXd G(n_assets, n_assets);//=mu*sigma
            VectorXd g0(n_assets);//-pbar
            MatrixXd CI = -1.0*MatrixXd::Identity(n_assets, n_assets);
            VectorXd ci0 = VectorXd::Zero(n_assets);

            MatrixXd CE = MatrixXd::Identity(n_assets, n_assets);
            VectorXd ce0(1);
            ce0 << -1.0; //weights sum to one constraint

            for (int p=0; p<n_portfolios; ++p){
                double mu = pow(10, 5.0*p/n_portfolios - 1.0); //log spacing of returns to query
                G = covariance * mu;
                g0 = -returns;

                VectorXd local_weights(n_assets);
                std::cout<<solve_quadprog(G, g0, CE, ce0, CI, ci0, local_weights)<<std::endl;
                weights.row(p) = local_weights;

                profile.push_back(
                    std::make_pair(
                        returns.dot(local_weights), //profit
                        sqrt(local_weights.dot(covariance*local_weights)) //risk
                    )
                );
            }
        }
        return weights;
    }

}