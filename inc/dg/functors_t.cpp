#include <iostream>


#include "functors.h"

#include "catch2/catch_all.hpp"

//TODO Add tests for other functors

TEST_CASE("Basic Functors")
{
    std::vector<double> xs( { 1.,2.,3.,4.});
    SECTION("Random")
    {
        dg::RandomNumbers<float> rand( 0,1);
        for( unsigned u=0; u<xs.size(); u++)
        {
            // Manually check that random numbers look random
            //std::cout << rand(xs[u], xs[u], xs[u], xs[u])<<"\n";

            CHECK( 0<=rand( xs[u], xs[u], xs[u]));
            CHECK( rand( ) <1);
        }
        // Test that rand can be called on device
        thrust::device_vector<float> xd( {1,2,3,4});
        thrust::transform( xd.begin(), xd.end(), xd.begin(), rand);
        for( unsigned u=0; u<xd.size(); u++)
        {
            CHECK( 0<=xd[u]);
            CHECK( xd[u] <1);
        }
    }
}

TEST_CASE( "Horner")
{
    // 10-the Legendre Polynomial
    std::vector<double> c10 = {-63, 0., 3465, 0, -30030, 0, 90090, 0, -109395, 0, 46189};
    dg::blas1::scal( c10, 1.0/256.);
    // 7-the Legendre Polynomial
    std::vector<double> c7 = {0., -35, 0, 315, 0, -693, 0, 429};
    dg::blas1::scal( c7, 1.0/16.);
    const double leg10 = -0.122124997387109375; // at x = 0.1
    const double dxleg10 = 2.258587288632812;  // at x = 0.1
    const double leg7 = -0.2935168; // at x = 0.2
    const double dxleg7 = -0.159488; // at x = 0.2
    SECTION("Horner1d")
    {
        // Test Legendre polynomial
        dg::Horner1d legendre( c10);
        double leg10_num = legendre(0.1);
        CHECK( fabs ( leg10_num - leg10) < 1e-12);
        auto dx = legendre.dx();
        CHECK( fabs( dx(0.1) - dxleg10) <1e-12);
    }
    SECTION("Horner2d")
    {
        // Multiply L10[x]*L7[y]
        std::vector<double> c107( c10.size()*c7.size());
        for( unsigned i=0; i<11; i++)
            for( unsigned j=0; j<8; j++)
                c107[i*8+j] = c10[i]*c7[j];
        dg::Horner2d horner( c107, c10.size(), c7.size());

        CHECK( fabs( horner( 0.1,0.2) - leg10*leg7) < 1e-12);
        CHECK( fabs( horner.dx()( 0.1,0.2) - dxleg10*leg7) < 1e-12);
        CHECK( fabs( horner.dy()( 0.1,0.2) - leg10*dxleg7) < 1e-12);
    }

}
