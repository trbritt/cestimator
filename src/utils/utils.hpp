/*
 *  author: Tristan Britt
 *  email: hello@tbritt.xyz
 *
 *  file: utils.hpp
 *  license: gpl_v3
 *
 *  this is the header file for all utilities of the project,
 *  namely anything not theoretically driving the main results,
 *  but necessary to get there
 *
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <chrono>
#ifdef __VISUALIZER
#include <matplot/matplot.h>
#endif
using namespace Eigen;

namespace Cestimator {
namespace Utils {
class Timer;

MatrixXd covariance(MatrixXd& x);

VectorXd mean(MatrixXd& x);

VectorXd outlier_cutoff(VectorXd& d, double d0);

void banner(std::string fname);

void goodbye(std::vector <Timer> timers);
void goodbye();

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
std::vector <double> matrix2vector(MatrixXd& x);

        #ifdef __VISUALIZER
class Visualizer {
public:
   Visualizer() {
      matplot::hold(matplot::on);
   };
   int scatter3d(MatrixXd& arr);
   int scatter2d(MatrixXd& arr, int dim1, int dim2);
   int ellipse(const VectorXd& mu, const MatrixXd& sigma, int dim1, int dim2);
};
        #endif
}
}

class Cestimator::Utils::Timer {
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
template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector <T> vec) {
   os << "{";
   if (vec.size() != 0) {
      std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator <T>(os, " "));
      os << vec.back();
   }
   os << "}";
   return os;
}
