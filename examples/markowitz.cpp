/*
 *  author: Tristan Britt
 *  email: hello@tbritt.xyz
 *
 *  file: markowitz.cpp
 *  license: gpl_v3
 *
 */
#include "cestimator.hpp"

std::pair <VectorXd, MatrixXd> market_projection(const VectorXd& exp_yr_comp_rets, const MatrixXd& cov_yr_comp_rets, const VectorXd& starting_prices, const double& horizon) {
   VectorXd mu    = exp_yr_comp_rets * horizon;
   MatrixXd sigma = cov_yr_comp_rets * horizon;

   VectorXd M = (mu.array() + 0.5 * sigma.diagonal().array()).exp();

   VectorXd exp_prices = starting_prices.asDiagonal() * M;

   MatrixXd S = (M * M.transpose()).array() * (sigma.array().exp() - MatrixXd::Ones(sigma.rows(), sigma.cols()).array()).array();

   MatrixXd cov_prices = starting_prices.asDiagonal() * S * starting_prices.asDiagonal();
   return std::make_pair(exp_prices, cov_prices);
}

int main(int argc, char *argv[]) {
   if (argc != 2) {
      std::cerr << "Incorrect number of command line arguments passed to programme. Ending" << std::endl;
      return 1;
   }

   Eigen::initParallel();
   Eigen::setNbThreads(6);

   std::string fname = std::string(argv[1]);

   Cestimator::Utils::banner(fname);
   MatrixXd data = load_csv <MatrixXd>(fname);
   MatrixXd Y    = data(all, all);
   Y.transposeInPlace();

   MatrixXd X = (Y(all, seq(1, last)).array().log()).matrix() - (Y(all, seq(0, last - 1)).array().log()).matrix(); //weekly compounded returns

   MatrixXd expected_yr_compounded_returns   = X.rowwise().mean() * 52;
   MatrixXd covariance_yr_compounded_returns = Cestimator::Utils::covariance(X) * 52;

   VectorXd current_prices = Y(all, last);

   // For the sake of this example, we are assuming that the
   // markets are normally distributed, which makes their projection
   // to an arbitrary horizon straightforward, following the "sqrt root" rule
   double horizon = 1.0;
   std::pair <VectorXd, MatrixXd> projected_market = market_projection(expected_yr_compounded_returns, covariance_yr_compounded_returns, current_prices, horizon);

   VectorXd exp_prices = std::get <0>(projected_market);
   MatrixXd cov_prices = std::get <1>(projected_market);

   Cestimator::Optimizer random(Cestimator::Optimizer::Classifier::Qp);

   random.generate_frontier(exp_prices, cov_prices, 100);

   for (const auto& [first, sec, third] : random.profile) {
      std::cout << "Weights [%]:" << 100 * first.transpose() << std::endl;
      std::cout << "\t Return:" << sec << ",\tRisk:" << third << std::endl;
   }

   Cestimator::Utils::goodbye();
   return 0;
}
