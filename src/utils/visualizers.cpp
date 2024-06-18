/*
 *  author: Tristan Britt
 *  email: hello@tbritt.xyz
 *
 *  file: visualizers.cpp
 *  license: gpl_v3
 *
 *  this is the implementation for all visualization routines of the project,
 *
 */
#ifdef __VISUALIZER
#include "utils.hpp"

using namespace Eigen;
using namespace matplot;
namespace Cestimator {
namespace Utils {
std::vector<double> vec2stdvec(const VectorXd& a){
   std::vector<double> vec(a.data(), a.data() + a.size()); 
   return vec;
};

int Visualizer::scatter3d(MatrixXd& arr) {
   auto l = plot3(vec2stdvec(arr.row(0)), vec2stdvec(arr.row(1)), vec2stdvec(arr.row(2)), "-ob");
   l->marker_size(5);
   l->marker_face_color("#FFFFFF");
   xlabel("x_1(t)");
   ylabel("x_2(t)");
   zlabel("x_3(t)");
   return 0;
}

int Visualizer::scatter2d(MatrixXd& arr, int dim1, int dim2) {
   auto l = plot(vec2stdvec(arr.row(dim1)), vec2stdvec(arr.row(dim2)), "-ob");
   l->marker_size(5);
   xlabel("x_1(t)");
   ylabel("x_2(t)");
   return 0;
}

int Visualizer::ellipse(const VectorXd& mu, const MatrixXd& sigma, int dim1, int dim2) {
   EigenSolver <MatrixXd> sigma_solver(sigma({ dim1, dim2 }, { dim1, dim2 }));
   VectorXd evals   = sigma_solver.eigenvalues().real();
   MatrixXd devals  = evals.array().sqrt().matrix().asDiagonal();
   MatrixXd evecs   = sigma_solver.eigenvectors().real();
   int      n_angs  = 500;
   VectorXd ellipse_x = VectorXd::Zero(n_angs);
   VectorXd ellipse_y = VectorXd::Zero(n_angs);

   for (int i = 0; i < n_angs; ++i) {
      VectorXd y(2);
      y << cos(2 * M_PI * i / n_angs), sin(2 * M_PI * i / n_angs);
      ellipse_x(i) = (2 * evecs * devals * y)(0) + mu(0);
      ellipse_y(i) = (2 * evecs * devals * y)(1) + mu(1);
   }
   auto l = plot(vec2stdvec(ellipse_x), vec2stdvec(ellipse_y));
   l->line_width(10);
   // mglData x(n_angs), y(n_angs);
   // for (int i = 0; i < n_angs; ++i) {
   //    x.a[i] = ellipse(i, 0) + mu(0);
   //    y.a[i] = ellipse(i, 1) + mu(1);
   // }
   // gr->Plot(x, y);

   return 0;
}
}
}

#endif
