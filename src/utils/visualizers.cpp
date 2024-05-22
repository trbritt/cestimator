/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: visualizers.cpp
    license: gpl_v3

    this is the implementation for all visualization routines of the project,

*/
#ifdef __VISUALIZER
#include "utils.hpp"

using namespace Eigen;

namespace Cestimator{
    
    namespace Utils {

        int Visualizer::scatter3d(MatrixXd& arr){
            int i, n=arr.cols();
            mglData x(n),y(n),z(n),c(n);
            for ( i = 0; i<n ;++i){
                x.a[i] = arr(0,i);
                y.a[i] = arr(1,i);
                z.a[i] = arr(2,i);
                c.a[i] = 0.;
            }
            gr->SetSize(2000,2000);
            gr->Alpha(true);
            gr->SubPlot(1,1,0); gr->Title("Scatter 3D"); gr->Rotate(50,60);
            gr->Box();  gr->Dots(x,y,z,c);
            return 0;
        }

        int Visualizer::scatter2d(MatrixXd& arr, int dim1, int dim2){
            int i, n=arr.cols();
            mglData x(n),y(n),r(n);
            for ( i = 0; i<n ;++i){
                x.a[i] = arr(dim1,i);
                y.a[i] = arr(dim2,i);
                r.a[i] = 0.5;
            }
            gr->SetSize(2000,2000);
            gr->Alpha(true);
            gr->SubPlot(1,1,0); gr->Title("Scatter 2D"); 
            gr->Box();  
            gr->Mark(x,y,r,"ko");
            
            return 0;
        }

        int Visualizer::ellipse(const VectorXd& mu, const MatrixXd& sigma, int dim1, int dim2){
            EigenSolver<MatrixXd> sigma_solver(sigma);
            VectorXd evals = sigma_solver.eigenvalues().real();
            MatrixXd devals = evals({dim1,dim2}).array().sqrt().matrix().asDiagonal();
            MatrixXd evecs = sigma_solver.eigenvectors().real();
            int n_angs = 500;
            MatrixXd ellipse = MatrixXd::Zero(n_angs,2);
            for (int i=0; i<n_angs; ++i){
                VectorXd y(2);
                y << cos(2*M_PI*i/n_angs), sin(2*M_PI*i/n_angs);
                ellipse.row(i) = 2*evecs({0,2}, {0,2}) * devals * y;
            }
            
            mglData x(n_angs),y(n_angs);
            for (int i = 0; i<n_angs ;++i){
                x.a[i] = ellipse(i,0)+mu(0);
                y.a[i] = ellipse(i,1)+mu(1);
            }
            gr->Plot(x, y, "r");

            return 0;
        }
    }
}

#endif
