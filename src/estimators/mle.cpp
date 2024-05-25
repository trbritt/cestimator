/*
 *  author: Tristan Britt
 *  email: hello@tbritt.xyz
 *
 *  file: mle.cpp
 *  license: gpl_v3
 *
 *  this is the main code for the maximum likelihood estimator
 *
 */
#include "estimators.hpp"
#include <random>

template <typename T>
static inline double lerp(T v0, T v1, T t) {
   return (1 - t) * v0 + t * v1;
}

template <typename T>
static inline std::vector <T> quantile(const std::vector <T>& inData, const std::vector <T>& probs) {
   if (inData.empty()) {
      return std::vector <T>();
   }

   if (1 == inData.size()) {
      return std::vector <T>(1, inData[0]);
   }

   std::vector <T> data = inData;
   std::sort(data.begin(), data.end());
   std::vector <T> quantiles;

   for (size_t i = 0; i < probs.size(); ++i) {
      T poi = lerp <T>(-0.5, data.size() - 0.5, probs[i]);

      size_t left  = std::max(int64_t(std::floor(poi)), int64_t(0));
      size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));

      T datLeft  = data.at(left);
      T datRight = data.at(right);

      T quantile = lerp <T>(datLeft, datRight, poi - left);

      quantiles.push_back(quantile);
   }

   return quantiles;
}

int Cestimator::maximum_likelihood::run() noexcept{
   double error = 1e6;


   VectorXd w     = VectorXd::Ones(T);
   VectorXd Zeros = VectorXd::Zero(N);

   mu = Zeros;

   sigma = MatrixXd::Zero(N, N);

   std::vector <double> flattened = Cestimator::Utils::matrix2vector(data);
   std::vector <double> quant     = quantile <double>(flattened, { 0.75, 0.25 });

   const double tolerance = abs(0.01 * (quant[1] - quant[0]));

   const int nus[6] = { 1, 2, 4, 7, 12, 20 };

   std::vector <double>             LL;
   std::vector <Cestimator::Result> res;
   VectorXd mu_old    = VectorXd::Zero(N);
   MatrixXd sigma_old = MatrixXd::Zero(N, N);
   MatrixXd inv_sigma = MatrixXd::Identity(N, N);
   MatrixXd W(N, T);
   MatrixXd x_c(N, T);
   VectorXd ma2(T);
   MatrixXd tmp(N, N);

   for (int idn = 0; idn < 6; ++idn) {
      double nu = static_cast <double>(nus[idn]);
      while (error > tolerance) {
         mu_old    = mu;
         sigma_old = sigma;

         mu        = (data * (w.array().square()).matrix()).rowwise().sum() / w.sum();
         W         = (data.colwise() - mu) * w.asDiagonal();
         sigma     = W * W.transpose() / (w.array().square().sum());
         x_c       = data.colwise() - mu;
         inv_sigma = sigma.llt().solve(MatrixXd::Identity(N, N));
         ma2       = (x_c.transpose() * inv_sigma * x_c).rowwise().sum();
         w         = (nu + N) / (nu + ma2.array());
         tmp       = (sigma - sigma_old).cwiseProduct(sigma - sigma_old) / N;
         tmp      += (mu - mu_old) * (mu - mu_old).transpose() / N;
         error     = tmp.trace();
      }
      // now that we've achieved the mu and sigma for
      // this given degree of freedom, compute its LL
      double norm = -N / 2 * log(nu * M_PI) + std::lgamma((nu + N) / 2) - std::lgamma(nu / 2) - 0.5 * log(sigma.determinant());

      double ll = 0;
      for (int t = 0; t < T; ++t) {
         MatrixXd centered = data.col(t) - mu;
         double   ma2      = (centered.transpose() * inv_sigma * centered)(0);
         ll += norm - (nu + N) / 2 * log(1 + ma2 / nu);
      }
      res.push_back(std::make_tuple(mu, sigma));
      LL.push_back(ll);
   }
   std::vector <double>::iterator result = std::max_element(LL.begin(), LL.end());
   int argmaxVal             = std::distance(LL.begin(), result);
   int nu                    = nus[argmaxVal];
   Cestimator::Result retval = res[argmaxVal];
   mu     = std::get <0>(retval);
   sigma  = std::get <1>(retval);
   sigma *= nu / (nu - 2);
   return 0;
}

int Cestimator::GMM::run() noexcept{
   const int iterations = 1000;

   int N = data.rows();  //number of dimensions
   int T = data.cols();  //number of observations

   std::vector <double> LL;
   double previous_LL = -std::numeric_limits <double>::infinity();

   mu = MatrixXd::Zero(n_features, N);
   VectorXd pi = VectorXd::Ones(n_features) / n_features;

   sigma = std::vector <MatrixXd>(n_features, MatrixXd::Identity(N, N));

   regularization = 1e-6 * MatrixXd::Identity(N, N);

   std::random_device rd;
   std::mt19937       gen(rd());
   std::uniform_int_distribution <> dis(0, T - 1);

   // Initialize mu
   for (int c = 0; c < n_features; ++c) {
      mu.row(c) = data.col(dis(gen));
   }

   // Initialize covariances
   for (int c = 0; c < n_features; ++c) {
      sigma[c] = 5 * MatrixXd::Identity(N, N);
   }

   for (int iter = 0; iter < iterations; ++iter) {
      // E Step
      MatrixXd r_tc = MatrixXd::Zero(T, n_features);
      for (int c = 0; c < n_features; ++c) {
         MatrixXd sigma_local = sigma[c] + regularization;
         for (int i = 0; i < T; ++i) {
            VectorXd diff     = data.col(i) - mu.row(c).transpose();
            double   exponent = -0.5 * diff.transpose() * sigma_local.llt().solve(MatrixXd::Identity(N, N)) * diff;
            double   denom    = pow(2 * M_PI, N / 2.0) * sqrt(sigma_local.determinant());
            r_tc(i, c) = pi(c) * exp(exponent) / denom;
         }
      }
      for (int t = 0; t < T; ++t) {
         r_tc.row(t) /= r_tc.row(t).sum();
      }
      // M Step
      for (int c = 0; c < n_features; ++c) {
         double m_c = r_tc.col(c).sum();
         mu.row(c) = (data * r_tc.col(c)).transpose() / m_c;

         MatrixXd covariance = MatrixXd::Zero(N, N);
         for (int t = 0; t < T; ++t) {
            VectorXd diff = data.col(t) - mu.row(c).transpose();
            covariance += r_tc(t, c) * diff * diff.transpose();
         }
         sigma[c] = covariance / m_c + regularization;
         pi(c)    = m_c / T;
      }

      // LL termination criterion
      double LL = 0;
      for (int t = 0; t < T; ++t) {
         double sum = 0.0;
         for (int c = 0; c < n_features; ++c) {
            MatrixXd sigma_local = sigma[c] + regularization;
            VectorXd diff        = data.col(t) - mu.row(c).transpose();
            double   exponent    = -0.5 * diff.transpose() * sigma_local.llt().solve(MatrixXd::Identity(N, N)) * diff;
            double   denom       = pow(2 * M_PI, data.rows() / 2.0) * sqrt(sigma_local.determinant());
            sum += pi(c) * exp(exponent) / denom;
         }
         LL += log(sum);
      }
      if (abs(LL - previous_LL) < 1e-4) {
         // std::cout << "Converged at iteration " << iter << std::endl;
         n_term = iter;
         break;
      }
      previous_LL = LL;
   }
   return 0;
}

VectorXd Cestimator::GMM::predict(VectorXd& arr) {
   VectorXd pred = VectorXd::Zero(n_features);
   for (int c = 0; c < n_features; ++c) {
      MatrixXd sigma_local = sigma[c] + regularization;
      VectorXd diff        = arr - mu.row(c).transpose();
      double   exponent    = -0.5 * diff.transpose() * sigma_local.llt().solve(MatrixXd::Identity(N, N)) * diff;
      double   denom       = pow(2 * M_PI, arr.size() / 2.0) * sqrt(sigma_local.determinant());
      pred(c) = exp(exponent) / denom;
   }
   pred /= pred.sum();
   return pred;
}

void Cestimator::GMM::print() {
   size_t      num_spaces = std::max(0, static_cast <int>(15 - name.length()));
   std::string padding(num_spaces, ' ');

   int ndims = data.rows();
   std::cout << Cestimator::Utils::colors::OKCYAN << name + padding << "\t\t" << "μ" << "\t\t" << "Σ" << Cestimator::Utils::colors::ENDC << std::endl;
   if (n_term == -1) {
      std::cout << Cestimator::Utils::colors::FAIL << "***" << name << " did not converge***" << Cestimator::Utils::colors::ENDC << std::endl;
   }else {
      std::cout << Cestimator::Utils::colors::WARNING << "***" << name << " converged after " << n_term << " iterations ***" << Cestimator::Utils::colors::ENDC << std::endl;
   }

   for (int c = 0; c < n_features; ++c) {
      for (int i = 0; i < ndims; ++i) {
         std::cout << "\t\t" << mu.row(c)(i) << "\t\t";
         for (int j = 0; j < ndims; ++j) {
            std::cout << sigma[c](i, j) << " ";
         }
         std::cout << std::endl;
      }
      std::cout << "\n" << std::endl;
   }
}
