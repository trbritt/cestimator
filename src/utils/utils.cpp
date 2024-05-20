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

namespace Cestimator{
    
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

        VectorXd mean(MatrixXd x){
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

        #ifdef __VISUALIZER
        using namespace Eigen;

        int scatter3d(mglGraph *gr, MatrixXd arr){
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

        int scatter2d(mglGraph *gr, MatrixXd arr, int dim1, int dim2){
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
        #endif
    }
}