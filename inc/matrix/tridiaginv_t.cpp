#include <iostream>
#include <iomanip>

#include "dg/algorithm.h"
#include "tridiaginv.h"

#include "catch2/catch_all.hpp"

TEST_CASE( "tridiaginv")
{
    std::vector<double> a = {1.98242,      4.45423,     5.31867,    7.48144,  7.11534};
    std::vector<double> b = {-0.00710891, -0.054661, -0.0554193, -0.0172191, -0.297645};
    std::vector<double> c = {0., -1.98242,     -4.44712,   -5.26401,   -7.42602, -7.09812};
    unsigned size = a.size();

    dg::TriDiagonal<thrust::host_vector<double>> T(size);
    dg::SquareMatrix<double> Tinv, Tinv_sol, Tinv_error( size, 0.);
    dg::SquareMatrix<double> H(size);
    H(0,0) = 0.505249; H(0,1) = 0.000814795; H(0,2) =  8.4358e-6;  H(0,3) = 6.26392e-8;   H(0,4) =  1.51587e-10;
    H(1,0) = 0.227217; H(1,1) = 0.227217;    H(1,2) =  0.00235244; H(1,3) = 0.0000174678; H(1,4) =  4.22721e-8;
    H(2,0) = 0.19139;  H(2,1) = 0.19139;     H(2,2) =  0.19139;    H(2,3) = 0.00142115;   H(2,4) =  3.43918e-6;
    H(3,0) = 0.134988; H(3,1) = 0.134988;    H(3,2) =  0.134988;   H(3,3) = 0.134988;     H(3,4) =  0.000326671;
    H(4,0) = 0.140882; H(4,1) = 0.140882;    H(4,2) =  0.140882;   H(4,3) = 0.140882;     H(4,4) = 0.140882;
    Tinv_sol = H;

    for( unsigned i=0; i<size; i++)
    {
        T.O[i]   =  a[i];  // 0 diagonal
        T.M[i]   =  c[i];  // -1 diagonal
        T.P[i]   =  b[i];  // +1 diagonal
    }
    INFO( "T matrix");
    for( unsigned u=0; u<size; u++)
        INFO( T.M[u]<<" "<<T.O[u]<<" "<<T.P[u]);

    SECTION( "TridiagInvDF")
    {
        dg::mat::TridiagInvDF<double> invtridiag(size);
        invtridiag(T, Tinv);
        Tinv_error = Tinv_sol - Tinv;
        INFO( "Difference between Tinv and solution "<<Tinv_error);
        for( unsigned u=0; u<size*size; u++)
            CHECK( fabs( Tinv_error.data()[u]) <1e-6); // The accuracy of the analytical solution is not better...
    }
    SECTION( "TridiagInvD")
    {
        dg::mat::TridiagInvD<double> invtridiag(size);
        invtridiag(T, Tinv);
        Tinv_error = Tinv_sol - Tinv;
        INFO( "Difference between Tinv and solution "<<Tinv_error);
        for( unsigned u=0; u<size*size; u++)
            CHECK( fabs( Tinv_error.data()[u]) <1e-6); // The accuracy of the analytical solution is not better...
    }
    SECTION( "TridiagInvHMGTI")
    {
        dg::mat::TridiagInvHMGTI<double> invtridiag(size);
        invtridiag(T, Tinv);
        Tinv_error = Tinv_sol - Tinv;
        INFO( "Difference between Tinv and solution "<<Tinv_error);
        for( unsigned u=0; u<size*size; u++)
            CHECK( fabs( Tinv_error.data()[u]) <1e-6); // The accuracy of the analytical solution is not better...
    }
    SECTION( "Size 1")
    {
        // P[0] and M[0] lie outside the matrix and must be ignored
        dg::TriDiagonal<thrust::host_vector<double>> T1(1);
        T1.O[0] = 4., T1.P[0] = 2., T1.M[0] = 3.;
        dg::mat::TridiagInvDF<double>( 1)( T1, Tinv);
        CHECK( Tinv(0,0) == 0.25);
        dg::mat::TridiagInvD<double>( 1)( T1, Tinv);
        CHECK( Tinv(0,0) == 0.25);
        dg::mat::TridiagInvHMGTI<double>( 1)( T1, Tinv);
        CHECK( Tinv(0,0) == 0.25);
        CHECK( dg::mat::invert( T1)(0,0) == 0.25);
    }
    SECTION( "invert")
    {
        auto Tinv = dg::mat::invert(T);
        Tinv_error = Tinv_sol - Tinv;
        INFO( "Difference between Tinv and solution "<<Tinv_error);
        for( unsigned u=0; u<size*size; u++)
            CHECK( fabs( Tinv_error.data()[u]) <1e-6); // The accuracy of the analytical solution is not better...
    }
}

TEMPLATE_TEST_CASE( "Lapack", "", double, float)
{
    using real = TestType;
    const double tol = std::is_same_v<real, double> ? 1e-14 : 1e-6;
    // Mathematica: Eigensystem[N[ToeplitzMatrix[PadRight[{2,-1}, 5]],16]]
    std::vector<double> lambda = {3.732050807568877, 3.000000000000000, 2.000000000000000, 1.000000000000000, 0.2679491924311227};
    // Eigenvectors are analytical only up to a scalar!!
    std::vector<double> evecs = {
        -0.2886751345948129, 0.5000000000000000, -0.5773502691896258, 0.5000000000000000, -0.2886751345948129,
    -0.5000000000000000, 0.5000000000000000, 0.0, -0.5000000000000000, 0.5000000000000000,
    -0.5773502691896258, 0.0, 0.5773502691896258, 0.0, -0.5773502691896258,
    0.5000000000000000, 0.5000000000000000, 0.0, -0.5000000000000000, -0.5000000000000000,
    -0.2886751345948129, -0.5000000000000000, -0.5773502691896258, -0.5000000000000000, -0.2886751345948129};


    dg::TriDiagonal<thrust::host_vector<real>> T(5);
    for( unsigned u=0; u<5; u++)
    {
        T.O[u] = 2;
        T.P[u] = T.M[u] = -1;
    }
    SECTION( "stev N")
    {
        dg::TriDiagonal<std::vector<real>> TT(T);
        dg::mat::lapack::stev( 'N', TT.O, TT.P, TT.M, TT.M);
        for( unsigned u=0; u<5; u++)
        {
            INFO( "Eigenvalue "<<u<<" "<<TT.O[u]);
            CHECK( fabs( TT.O[u] - lambda[4-u] ) < tol);
        }
    }
    SECTION( "stev V")
    {
        dg::TriDiagonal<std::vector<real>> TT(T);
        dg::SquareMatrix<real> evs( 5);
        dg::SquareMatrix<double> sol( evecs.begin(), evecs.end());
        std::vector<real> work ( 2*5-2);
        // In LAPACK_COL_MAJOR the eigenvectors become the rows of evs!
        dg::mat::lapack::stev( 'V', TT.O, TT.P, evs.data(), work);
        for( unsigned u=0; u<5; u++)
        {
            INFO( "Eigenvalue "<<u<<" "<<TT.O[u]);
            CHECK( fabs( TT.O[u] - lambda[4-u] ) < tol);
        }
        for( unsigned u=0; u<5; u++)
        {
            // An Eigenvector is determined only up to a scalar factor
            // c = <evs_u, sol_u> / <sol_u, sol_u>; stev normalises them
            double dot = 0, norm = 0, norm_evs = 0;
            for( unsigned k=0; k<5; k++)
            {
                dot  += evs(u,k)*sol(4-u,k);
                norm += sol(4-u,k)*sol(4-u,k);
                norm_evs += evs(u,k)*evs(u,k);
            }
            double c = dot/norm;
            INFO( "Eigenvector "<<u<<" factor "<<c);
            CHECK( fabs( sqrt( norm_evs) - 1.) < tol);
            for( unsigned k=0; k<5; k++)
            {
                INFO( "Eigenevector "<<u<<" "<<k<<" "<<evs(u,k));
                CHECK ( fabs( evs(u,k) - c*sol( 4-u,k)) < tol);
            }
        }
    }
    SECTION( "Compute extreme EV")
    {
        auto EVs = dg::mat::compute_extreme_EV( T);
        CHECK( fabs(EVs[0] - lambda[4]) < tol);
        CHECK( fabs(EVs[1] - lambda[0]) < tol);
    }
    SECTION( "sygv")
    {
        // A x = lambda B x with B = 2 I has Eigenvalues lambda/2
        std::vector<real> A( 25, 0), B( 25, 0), W( 5), work( 3*5-1);
        for( unsigned u=0; u<5; u++)
        {
            A[u*5+u] = 2, B[u*5+u] = 2;
            if( u < 4)
                A[u*5+u+1] = A[(u+1)*5+u] = -1;
        }
        dg::mat::lapack::sygv( 1, 'N', 'U', 5, A, 5, B, 5, W, work);
        for( unsigned u=0; u<5; u++)
        {
            INFO( "Eigenvalue "<<u<<" "<<W[u]);
            CHECK( fabs( W[u] - lambda[4-u]/2. ) < tol);
        }
    }
}
