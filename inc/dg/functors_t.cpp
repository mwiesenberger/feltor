#include <iostream>
#include <cmath>
#include <array>


#include "functors.h"
#include "blas2.h"
#include "topology/weights.h"

#include "catch2/catch_all.hpp"

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

TEST_CASE( "Product multiplies from left to right")
{
    // For these values the order of multiplication changes the result
    double x = 0.1, y = 0.2, z = 0.3, w = 0.4;
    REQUIRE( (x*y)*z != x*(y*z));
    CHECK( dg::Product()( x, y, z) == (x*y)*z);
    CHECK( dg::Product()( x, y, z, w) == ((x*y)*z)*w);
    CHECK( dg::Product()( x, y) == x*y);
    CHECK( dg::Product()( x) == x);
    // The result type follows the usual promotion of x*y
    CHECK( std::is_same_v<decltype( dg::Product()( 0.1f, 2, 0.3)), double>);
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
    SECTION( "Horner1d is correct when the point repeats")
    {
        // f(x) = 1 + 2x + 3x^2
        dg::Horner1d p( {1., 2., 3.});
        auto f = []( double x){ return 1. + x*( 2. + 3.*x);};
        for( double x : { 0.5, 0.5, -0.7, 0.5, 0., -0.7})
        {
            INFO( "x "<<x);
            CHECK( fabs( p( x) - f( x)) < 1e-15);
        }
        CHECK( dg::Horner1d()( 0.3) == 1.);
        CHECK( dg::Horner1d( std::vector<double>{})( 0.3) == 0.);
    }
    SECTION( "Horner2d is correct when the point repeats or changes in one coordinate")
    {
        // c[i*N+j] x^i y^j with M = 2, N = 3
        dg::Horner2d p( {1., 2., 3., 4., 5., 6.}, 2, 3);
        auto f = []( double x, double y){
            return 1. + 2.*y + 3.*y*y + x*( 4. + 5.*y + 6.*y*y);};
        // Same point, new x with same y, new y with same x, back to the start
        std::vector<std::array<double,2>> points = { {0.5, 0.3}, {0.5, 0.3},
            {-0.7, 0.3}, {-0.7, -1.1}, {0.5, -1.1}, {0.5, 0.3}, {0., 0.}};
        for( auto pt : points)
        {
            INFO( "x "<<pt[0]<<" y "<<pt[1]);
            CHECK( fabs( p( pt[0], pt[1]) - f( pt[0], pt[1])) < 1e-14);
        }
        dg::Horner2d one;
        CHECK( one( 0.3, 0.4) == 1.);
        CHECK( one( 0.3, -0.4) == 1.);
        CHECK( dg::Horner2d( {}, 0, 3)( 0.3, 0.4) == 0.);
        CHECK( dg::Horner2d( {}, 2, 0)( 0.3, 0.4) == 0.);
    }
}

TEST_CASE( "MOD returns values in [0,m)")
{
    dg::MOD<double> mod( 2.);
    CHECK( mod( 5.5) == 1.5);
    CHECK( mod( -0.5) == 1.5);
    CHECK( mod( 4.) == 0.);
    CHECK( mod( -4.) == 0.);
    CHECK( mod( 1.) == 1.);
    // Results that round to m must not be returned
    for( double x : { -1e-20, -1e-300, -4.-1e-17, -2.*(1.-1e-17), 0., 1.999})
    {
        INFO( "x "<<x<<" mod "<<mod( x));
        CHECK( 0. <= mod( x));
        CHECK( mod( x) < 2.);
    }
    dg::MOD<float> modf( 1.f);
    CHECK( 0.f <= modf( -1e-10f));
    CHECK( modf( -1e-10f) < 1.f);
}

TEST_CASE( "VanLeer limiter properties")
{
    dg::VanLeer vl;
    CHECK( vl( 1., -2.) == 0.);
    CHECK( vl( 0., 3.) == 0.);
    CHECK( vl( 3., 0.) == 0.);
    for( auto ab : std::vector<std::array<double,2>>{ {1., 2.}, {0.1, 5.},
            {3., 3.}, {-1., -4.}})
    {
        double a = ab[0], b = ab[1];
        INFO( "a "<<a<<" b "<<b);
        CHECK( vl( a, b) == vl( b, a));
        // Between the smaller and twice the smaller argument (in magnitude)
        double mn = std::min( fabs( a), fabs( b));
        CHECK( mn <= fabs( vl( a, b)));
        CHECK( fabs( vl( a, b)) <= 2.*mn);
        CHECK( vl( a, b)*a > 0);
    }
    CHECK( vl( 3., 3.) == 3.);
}

TEST_CASE( "PolynomialRectangle")
{
    double xl = -1., al = 0.5, xr = 2., ar = 0.3;
    dg::PolynomialRectangle rect( xl, al, xr, ar);
    CHECK( fabs( rect( xl) - 0.5) < 1e-14);
    CHECK( fabs( rect( xr) - 0.5) < 1e-14);
    CHECK( rect( xl - 1.1*al) == 0.);
    CHECK( rect( xr + 1.1*ar) == 0.);
    CHECK( rect( 0.5*(xl+xr)) == 1.);
    CHECK( rect( xl + 1.1*al) == 1.);
    CHECK( rect( xr - 1.1*ar) == 1.);
    // Monotone in both transition regions
    for( unsigned i=0; i<20; i++)
    {
        double x = xl - al + 2.*al*i/20.;
        CHECK( rect( x) <= rect( x + al/20.));
        x = xr - ar + 2.*ar*i/20.;
        CHECK( rect( x) >= rect( x + ar/20.));
    }
    // With equal widths the rectangle is symmetric around its centre
    dg::PolynomialRectangle sym( xl, al, xr, al);
    for( double x : { -1.3, -0.8, 1.6, 2.4})
    {
        INFO( "x "<<x);
        CHECK( fabs( sym( x) - sym( xl+xr-x)) < 1e-14);
    }
}

TEST_CASE( "Lamb dipole")
{
    const double x0 = 0.2, y0 = -0.1, R = 1., U = 0.5;
    dg::Lamb lamb( x0, y0, R, U);
    // The vorticity vanishes at the edge of the dipole
    for( double theta : { 0., 0.4, 2.})
    {
        INFO( "theta "<<theta);
        double x = x0 + R*cos( theta)*(1.-1e-14), y = y0 + R*sin( theta)*(1.-1e-14);
        CHECK( fabs( lamb( x, y)) < 1e-12);
        CHECK( lamb( x0 + 1.1*R*cos( theta), y0 + 1.1*R*sin( theta)) == 0.);
    }
    // enstrophy = 1/2 int omega^2
    dg::Grid2d grid( -2., 2., -2., 2., 3, 100, 100);
    auto omega = dg::evaluate( lamb, grid);
    double enstrophy = 0.5*dg::blas2::dot( omega, dg::create::weights( grid),
        omega);
    INFO( "Numerical "<<enstrophy<<" analytical "<<lamb.enstrophy());
    CHECK( fabs( enstrophy - lamb.enstrophy()) < 1e-6*lamb.enstrophy());
}

TEST_CASE( "Vortex")
{
    const double R = 1.5, ud = 0.7;
    unsigned state = GENERATE( 0, 1, 2);
    INFO( "State "<<state);
    dg::Vortex v( 0., 0., state, R, ud);
    // f(r,0) is continuous and continuously differentiable at r = R
    auto f = [&]( double r){ return v( r, 0.);};
    double h = 1e-4;
    double fin = f( R*(1.-1e-14)), fout = f( R*(1.+1e-14));
    CHECK( fabs( fin - fout) < 1e-10*fabs( fin));
    // one-sided second order differences; the constants g_i, beta_i have 10
    // digits and g_i lies close to a zero of J_1, so the slopes agree to 1e-5
    double din  = ( 3.*f( R) - 4.*f( R-h) + f( R-2.*h))/2./h;
    double dout = (-3.*f( R) + 4.*f( R+h) - f( R+2.*h))/2./h;
    INFO( "Derivative inside "<<din<<" outside "<<dout);
    CHECK( fabs( din - dout) < 1e-5*fabs( din));
    // The angle dependence is cos(theta)
    for( double theta : { 0.3, 2.})
    {
        for( double r : { 0.5*R, 2.*R})
        {
            INFO( "theta "<<theta<<" r "<<r);
            CHECK( fabs( v( r*cos( theta), r*sin( theta)) - f( r)*cos( theta))
                    < 1e-14*fabs( f( r)));
        }
    }
#ifdef __cpp_lib_math_special_functions
    // Outside the vortex f(r,0)/f(R,0) = K_1( beta r/R) / K_1( beta)
    // with the documented beta_i
    double beta = state == 0 ? 0.03827327723 : 0.07071067810;
    for( double r : { 1.2*R, 3.*R, 50.*R, 200.*R})
    {
        INFO( "r "<<r);
        double ratio = std::cyl_bessel_k( 1., beta*r/R)/std::cyl_bessel_k( 1.,
                beta);
        CHECK( fabs( f( r)/f( R) - ratio) < 1e-7*ratio);
    }
#endif
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
