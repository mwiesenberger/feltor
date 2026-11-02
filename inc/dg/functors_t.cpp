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

TEST_CASE( "Derivatives of functors")
{
    // Central differences have an error of O(h^2)
    const double h = 1e-5, eps = 1e-7;
    auto diff = [h]( auto f, double x){ return (f(x+h)-f(x-h))/2./h;};
    SECTION( "Cauchy")
    {
        dg::Cauchy c( 0.1, -0.2, 1.3, 0.7, 2.);
        double x = 0.4, y = 0.1;
        double dx = diff( [&]( double t){ return c( t, y);}, x);
        double dy = diff( [&]( double t){ return c( x, t);}, y);
        double dxx = diff( [&]( double t){ return c.dx( t, y);}, x);
        double dyy = diff( [&]( double t){ return c.dy( x, t);}, y);
        double dxy = diff( [&]( double t){ return c.dx( x, t);}, y);
        CHECK( fabs( c.dx( x, y) - dx) < eps*fabs( dx));
        CHECK( fabs( c.dy( x, y) - dy) < eps*fabs( dy));
        CHECK( fabs( c.dxx( x, y) - dxx) < eps*fabs( dxx));
        CHECK( fabs( c.dyy( x, y) - dyy) < eps*fabs( dyy));
        CHECK( fabs( c.dxy( x, y) - dxy) < eps*fabs( dxy));
    }
    SECTION( "PolynomialHeaviside")
    {
        // IPolynomialHeaviside' = PolynomialHeaviside and
        // PolynomialHeaviside' = DPolynomialHeaviside
        int sign = GENERATE( +1, -1);
        double xb = 0.3, a = 0.5;
        dg::PolynomialHeaviside H( xb, a, sign);
        dg::IPolynomialHeaviside IH( xb, a, sign);
        dg::DPolynomialHeaviside DH( xb, a, sign);
        INFO( "Sign "<<sign);
        CHECK( fabs( H( xb) - 0.5) < 1e-15);
        CHECK( H( xb - 1.1*a) == ( sign == 1 ? 0. : 1.));
        CHECK( H( xb + 1.1*a) == ( sign == 1 ? 1. : 0.));
        // The integration constant: IH is xb on the damped side and x on the other
        CHECK( fabs( IH( xb - 1.1*a) - ( sign == 1 ? xb : xb - 1.1*a)) < 1e-15);
        CHECK( fabs( IH( xb + 1.1*a) - ( sign == 1 ? xb + 1.1*a : xb)) < 1e-15);
        CHECK( DH( xb - 1.1*a) == 0.);
        CHECK( DH( xb + 1.1*a) == 0.);
        for( double x : { 0.0, 0.25, 0.5, 0.7})
        {
            INFO( "x "<<x);
            CHECK( fabs( diff( IH, x) - H( x)) < eps);
            CHECK( fabs( diff( H, x) - DH( x)) < eps);
        }
    }
}
