/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: mle.cpp
    license: gpl_v3

    this is the main code for the maximum likelihood estimator

*/
#include "estimators.hpp"

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