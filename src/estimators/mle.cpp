/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: mle.cpp
    license: gpl_v3

    this is the main code for the maximum likelihood estimator

*/
#include "estimators.hpp"
#include <random>

std::random_device rd;     // Only used once to initialise (seed) engine
std::mt19937 rng(rd());    // Random-number engine used (Mersenne-Twister in this case)

template<typename T>
static inline double lerp(T v0, T v1, T t)
{
    return (1 - t)*v0 + t*v1;
}

template<typename T>
static inline std::vector<T> quantile(const std::vector<T>& inData, const std::vector<T>& probs)
{
    if (inData.empty())
    {
        return std::vector<T>();
    }

    if (1 == inData.size())
    {
        return std::vector<T>(1, inData[0]);
    }

    std::vector<T> data = inData;
    std::sort(data.begin(), data.end());
    std::vector<T> quantiles;

    for (size_t i = 0; i < probs.size(); ++i)
    {
        T poi = lerp<T>(-0.5, data.size() - 0.5, probs[i]);

        size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
        size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));

        T datLeft = data.at(left);
        T datRight = data.at(right);

        T quantile = lerp<T>(datLeft, datRight, poi - left);

        quantiles.push_back(quantile);
    }

    return quantiles;
}

std::vector<double> matrix2vector(MatrixXd x){
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

int Cestimator::maximum_likelihood::run(){

    double error = 1e6;


    VectorXd w = VectorXd::Ones(T);
    VectorXd Zeros = VectorXd::Zero(N);

    mu = Zeros;

    sigma = MatrixXd::Zero(N,N);

    std::vector<double> flattened = matrix2vector(data);
    std::vector<double> quant = quantile<double>(flattened, {0.75, 0.25});

    const double tolerance = abs(0.01*(quant[1]-quant[0]));

    const int nus[6] = {1, 2, 4, 7, 12, 20};

    std::vector<double> LL;
    std::vector<Cestimator::Result> res;
    VectorXd mu_old;
    MatrixXd sigma_old;
    MatrixXd inv_sigma;
    MatrixXd W;
    MatrixXd x_c;
    VectorXd ma2;
    MatrixXd tmp;

    for (int idn=0; idn<6; ++idn){
        double nu = static_cast<double>(nus[idn]);
        while (error > tolerance) {
            mu_old = mu;
            sigma_old = sigma;

            W = w * VectorXd::Ones(N).transpose();
            mu = data.transpose().cwiseProduct(W).colwise().sum() / w.sum();
            x_c = data - mu * VectorXd::Ones(T).transpose();
            sigma = W.transpose().cwiseProduct(x_c) * x_c.transpose() / T;

            inv_sigma = sigma.inverse();
            ma2 = (x_c.transpose() * inv_sigma * x_c).rowwise().sum();
            for (size_t i=0; i<w.size();  ++i){
                w(i) = (nu+N) / (nu + ma2(i));
            }
            tmp = (sigma-sigma_old)*(sigma-sigma_old) / N;
            tmp += (mu-mu_old)*(mu-mu_old).transpose()/ N;
            error = tmp.trace();
        }
        // now that we've achieved the mu and sigma for
        // this given degree of freedom, compute its LL
        double norm = -N/2 * log(nu*M_PI) + std::lgamma((nu+N)/2) - std::lgamma(nu/2) - 0.5*log(sigma.determinant());

        double ll=0;
        for (int t=0; t<T; ++t){
            MatrixXd centered = data(all, t) - mu;
            double ma2 = (centered.transpose() * inv_sigma * centered)(0);
            ll += norm - (nu+N)/2 * log(1+ma2/nu);
        }
        res.push_back(std::make_tuple(mu, sigma));
        LL.push_back(ll);
    }
    std::vector<double>::iterator result = std::max_element(LL.begin(), LL.end());
    int argmaxVal = std::distance(LL.begin(), result);
    int nu = nus[argmaxVal];
    Cestimator::Result retval = res[argmaxVal];
    mu = std::get<0>(retval);
    sigma = std::get<1>(retval);
    sigma *= nu / (nu-2);
    return 0;
}

int Cestimator::GMM(MatrixXd X, int n_features){
    const int iterations = 5000;


    int N = X.rows(); //number of dimensions
    int T = X.cols(); //number of observations
    std::uniform_int_distribution<int> uni(
        X(0,all).minCoeff(),
        X(0,all).maxCoeff()
    ); // Guaranteed unbiased

    //because it is possible for a scatter matrix to become
    //non-invertible if we find the wrong local min, we add
    //regulatization to ensure smooth behaviour 
    MatrixXd regularization = 1e-6 * MatrixXd::Identity(N,N);

    //because we're in N-dims, and we have n_features number of
    //sources, mu must be an N x n_features matrix.
    std::vector<VectorXd> mu(n_features);
    std::fill(mu.begin(), mu.end(), VectorXd::Zero(N));
    //likewise, sigma must be n_features x N x N
    //fill the scatter matrix of each feature with 
    //an initial (uncorrelated) diagonal matrix
    std::vector<MatrixXd> sigma(n_features);
    std::fill(sigma.begin(), sigma.end(), MatrixXd::Zero(N, N));
    for (int c=0; c<n_features; ++c){
        for (int j=0; j<N; ++j){
            mu[c](j) = static_cast<double>(uni(rng));
            sigma[c](j, j) = 5.0;
        }
        // std::cout << mu[c] << std::endl;
        // std::cout << sigma[c] << std::endl;
    }

    //now, we define the mixing ratios of each source
    VectorXd pi = VectorXd::Ones(n_features) / n_features;
    double LL = 0;
    double LL_old = 0;
    for (int i=0; i<iterations; ++i){
        // std::cout << i << std::endl;
        MatrixXd r_ic = MatrixXd::Zero(n_features, T); //prob that this observations belongs to this cluster

        //first up is the 'E' step
        //start with iterating over n_features
        for (int c=0; c<n_features; ++c){
            VectorXd mu_c = mu[c];
            MatrixXd sigma_c = sigma[c]+regularization;

            MatrixXd inv_sigma_c = sigma[c].inverse();
            double num, denom=0;
            for (int t=0; t<T; ++t){
                VectorXd centered = X(all, t) - mu_c;
                num = exp(-0.5*
                    centered.transpose() * inv_sigma_c * centered
                ) * pi(c);
                num /= sqrt(pow(2*M_PI, N)*pow(sigma_c.determinant(), 0.5));
                r_ic(c, t) = num;
                denom += num;
            }
            r_ic(c, all) /= denom;
        }
        // now we begin with the M step
        //X is 2,500
        double m_c;
        for (int c=0; c<n_features; ++c){
            m_c = r_ic(c, all).sum();
            std::cout << m_c << std::endl;
            VectorXd mu_c = VectorXd::Zero(N);
            MatrixXd sigma_c = MatrixXd::Zero(N,N);

            for (int t=0; t<T; ++t){
                mu_c += r_ic(c,t) * X(all, t).transpose();
            }
            mu_c /= m_c;
            mu[c] = mu_c;

            for (int t=0; t<T; ++t){
                MatrixXd tmp = (X(all,t) - mu_c);
                sigma_c += r_ic(c,t) * tmp * tmp.transpose();
            }
            sigma[c] = sigma_c + regularization;

            pi(c) = m_c / r_ic.sum();

        }

        // now, compute LL for this iteration
        for (int c=0; c<n_features; ++c){
            VectorXd mu_c = mu[c];
            MatrixXd sigma_c = sigma[c];

            MatrixXd inv_sigma_c = sigma[c].inverse();
            double num_t=0, num_c=0;
            for (int t=0; t<T; ++t){
                VectorXd centered = X(all, t) - mu_c;
                num_t = exp(-0.5*
                    centered.transpose() * inv_sigma_c * centered
                ) * pi(c);
                num_t /= sqrt(pow(2*M_PI, N)*pow(sigma_c.determinant(), 0.5));
                num_c += num_t;
            }
            LL += num_c;
        }
        LL = log(LL);
        if (abs(LL-LL_old) < 1e-8) {
            std::cout << "Converged!" << std::endl;
            break;
        } else {
            LL_old = LL;
            std::cout << LL << std::endl;
        }
    }
    
    std::cout << mu << std::endl;
    // std::cout << sigma << std::endl;
    return 0;
};