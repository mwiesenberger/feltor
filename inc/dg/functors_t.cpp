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
    SECTION("Horner1d")
    {
        // Test Legendre polynomial
        std::vector<double> c = {-63, 0., 3465, 0, -30030, 0, 90090, 0, -109395, 0, 46189};
        dg::blas1::scal( c, 1.0/256.);
        dg::Horner1d legendre( c);
        const double leg10 = -0.122124997387109375;
        double leg10_num = legendre(0.1);
        CHECK( fabs ( leg10_num - leg10) < 1e-12);
    }

}

TEST_CASE( "Composition")
{
    SECTION( "Two functors")
    {
        // f(g(x,y)) = 2*(x+y)
        auto fg = dg::compose( []( double x){ return 2*x;},
                []( double x, double y){ return x+y;});
        CHECK( fg( 1., 2.) == 6.);
    }
    SECTION( "Three functors")
    {
        // f0(f1(f2(x,y))) = 2*(x*y)+1, the outermost functor is applied last
        auto f012 = dg::compose( []( double x){ return x+1;},
                []( double x){ return 2*x;},
                []( double x, double y){ return x*y;});
        CHECK( f012( 2., 3.) == 13.);
    }
}

TEST_CASE( "Histogram")
{
    SECTION( "Histogram")
    {
        // 4 bins [0,1), [1,2), [2,3), [3,4); values outside are discarded
        dg::Grid1d g1d( 0, 4, 1, 4);
        dg::Histogram<> hist( g1d, {0.5, 1.5, 1.6, 2.5, 2.6, 2.7, -0.5, 4.0, 5.});
        thrust::host_vector<double> centres = dg::evaluate( dg::cooX1d, g1d);
        std::vector<double> sol = { 1./3., 2./3., 1., 0.};
        for( unsigned u=0; u<4; u++)
        {
            INFO( "Bin "<<u);
            CHECK( hist( centres[u]) == sol[u]);
        }
        CHECK( hist( -0.5) == 0.);
        CHECK( hist( 4.0) == 0.);
        CHECK_THROWS_AS( dg::Histogram<>( dg::Grid1d( 0, 4, 2, 4), {0.5}),
            dg::Error);
    }
    SECTION( "Histogram2D")
    {
        // 2x2 bins; the pairs outside the domain are discarded
        dg::Grid2d g2d( 0, 2, 0, 2, 1, 2, 2);
        dg::Histogram2D<> hist( g2d, {0.5, 1.5, 1.5, 1.5, -1., 0.5},
                                     {0.5, 0.5, 1.5, 1.5, 0.5, 2.0});
        CHECK( hist( 0.5, 0.5) == 0.5);
        CHECK( hist( 1.5, 0.5) == 0.5);
        CHECK( hist( 0.5, 1.5) == 0.);
        CHECK( hist( 1.5, 1.5) == 1.);
        CHECK( hist( -1., 0.5) == 0.);
        CHECK( hist( 0.5, 2.0) == 0.);
        CHECK_THROWS_AS( dg::Histogram2D<>( dg::Grid2d( 0, 2, 0, 2, 2, 2, 2),
            {0.5}, {0.5}), dg::Error);
    }
}
