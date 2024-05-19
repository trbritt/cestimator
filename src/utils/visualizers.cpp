#include "utils.hpp"

#ifdef __VISUALIZER
using namespace Eigen;

int Cestimator::Utils::Visualizers::scatter3d(mglGraph *gr, MatrixXd arr){
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

int Cestimator::Utils::Visualizers::scatter2d(mglGraph *gr, MatrixXd arr, int dim1, int dim2){
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