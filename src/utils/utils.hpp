/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: utils.hpp
    license: gpl_v3

    this is the header file for all utilities of the project,
    namely anything not theoretically driving the main results,
    but necessary to get there

*/

#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <chrono>

namespace Utils {
    using namespace Eigen;

    class Timer;

    MatrixXd covariance(MatrixXd x);

    Vector2d mean(MatrixXd x);

    VectorXd outlier_cutoff(VectorXd d, double d0);

    void banner(std::string fname, int dim1, int dim2);

    void goodbye(std::vector<Timer> timers);

    std::string word_wrap(std::string text, unsigned per_line);

    class colors {
    public:
        static const std::string HEADER;
        static const std::string OKBLUE;
        static const std::string OKCYAN;
        static const std::string OKGREEN;
        static const std::string WARNING;
        static const std::string FAIL;
        static const std::string ENDC;
        static const std::string BOLD;
        static const std::string UNDERLINE;
    };
}

class Utils::Timer {
    public:
        typedef std::chrono::high_resolution_clock Clock;
        std::string name;
        void start();
        void stop();
        void print();
    private:
        Clock::time_point epoch;
        Clock::duration telapsed;
};

//pretty print vector
template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{";
    if(vec.size()!=0)
    {
        std::copy(vec.begin(), vec.end()-1, std::ostream_iterator<T>(os, " "));
        os<<vec.back();
    }
    os<<"}";
    return os;
}

