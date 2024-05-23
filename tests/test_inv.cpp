#include <gtest/gtest.h>

#define private public
#define protected public

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <chrono>
#include <iostream>

using namespace std::chrono;
using namespace Eigen;

typedef float T;
typedef high_resolution_clock Clock;

TEST(Cholesky, inverseTiming) {
    Matrix<T, Dynamic, Dynamic> L;
    Matrix<T, Dynamic, Dynamic> S;
    Matrix<T, Dynamic, Dynamic> Sinv_method1;
    Matrix<T, Dynamic, Dynamic> Sinv_method2;
    int Nmin = 2;
    int Nmax = 4;
    int N(Nmin);

    while (N <= Nmax) {

        L.resize(N, N);
        L.setRandom();
        S.resize(N, N);
        // create a random NxN SPD matrix
        S = L*L.transpose();
        std::cout << "\n";
        std::cout << "+++++++++++++++++++++++++ N = " << N << " +++++++++++++++++++++++++++++++++++++++" << std::endl;
        auto t1 = Clock::now();
        Sinv_method1.resize(N, N);
        Sinv_method1 = S.inverse();
        auto dt1 = Clock::now() - t1;
        std::cout << "Method 1 took " << duration_cast<microseconds>(dt1).count() << " usec" << std::endl;
        auto t2 = Clock::now();
        Sinv_method2.resize(N, N);
        Sinv_method2 = S.llt().solve(Matrix<T, Dynamic, Dynamic>::Identity(N, N));
        auto dt2 = Clock::now() - t2;
        std::cout << "Method 2 took " << duration_cast<microseconds>(dt2).count() << " usec" << std::endl;
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                EXPECT_NEAR( Sinv_method1(i, j), Sinv_method2(i, j), 1e-3 );
            }
        }

        N *= 2;
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "\n";
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}