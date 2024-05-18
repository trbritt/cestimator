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

    Vector2d mean(MatrixXd x){
        return x.rowwise().mean();
    }

    MatrixXd covariance(MatrixXd x){
        MatrixXd tmp = x.transpose();
        MatrixXd centered = tmp.rowwise() - tmp.colwise().mean();
        return (tmp.adjoint() * centered) / (double(tmp.rows() - 1));
    }

    VectorXd outlier_cutoff(VectorXd d, double d0){
        Matrix<int, -1, 1> index{d.size()};
        for (size_t i=0; i < d.size(); ++i){
            index(i) = d(i) <= d0 ? 1 : 0;
        }
        const double b_2 = 1.25;
        VectorXd omega{d.size()};
        for (size_t i = 0; i < d.size(); ++i){
            if (index(i)) {
                omega(i) = 1;
            } else {
                omega(i) = exp(-0.5*((d(i)-d0)*(d(i)-d0) / (b_2*b_2)));
            }
        }
        return omega;
    }

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
}