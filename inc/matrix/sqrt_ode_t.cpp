#include <iostream>

#include <boost/math/special_functions.hpp>
#include "sqrt_ode.h"

#include "catch2/catch_all.hpp"

using Container = dg::HVec;

// Multiplication with a diagonal d has analytic matrix functions
static const std::vector<double> diag_values = {0.5, 1., 2., 3.};
static const Container diag( diag_values.begin(), diag_values.end());
static void diagonal( const Container& y, Container& z)
{
    dg::blas1::pointwiseDot( diag, y, z);
}

TEST_CASE( "Bessel ODE")
{
    // z(1) = ( I_0(d) y(0), d I_1(d) y(0) ) for z(0) = ( y(0), 0)
    std::array<Container,2> z{ Container( diag.size(), 1.), Container( diag.size(), 0.)};
    std::array<Container,2> result( z);
    unsigned number = 0;
    auto ode = dg::mat::make_besselI0ode( diagonal);
    auto bessel = dg::mat::make_directODESolve( ode,
            std::string("Dormand-Prince-7-4-5"), 1e-10, 1e-12, number);
    bessel( z, result);
    for( unsigned u=0; u<diag.size(); u++)
    {
        INFO( "d "<<diag[u]);
        CHECK( fabs( result[0][u] - boost::math::cyl_bessel_i( 0, diag[u]))
                < 1e-8);
        CHECK( fabs( result[1][u] - diag[u]*boost::math::cyl_bessel_i( 1,
                diag[u])) < 1e-8);
    }
}

TEST_CASE( "InvSqrtODE")
{
    // y(1) = d^{-1/2} y(0)
    dg::mat::InvSqrtODE<Container> ode;
    ode.construct( diagonal, diag);
    ode.set_inverse_operator( [&t = ode.time()]( const Container& x, Container& y)
    {
        // ( t + (1-t) d)^{-1}
        dg::blas1::evaluate( y, dg::equals(), [&t]( double xx, double dd){
                return xx/( t + (1.-t)*dd);}, x, diag);
    });
    Container y0( diag.size(), 1.), y1( y0);
    unsigned number = 0;
    auto sqrtinv = dg::mat::make_directODESolve( ode,
            std::string("Dormand-Prince-7-4-5"), 1e-10, 1e-12, number);
    sqrtinv( y0, y1);
    for( unsigned u=0; u<diag.size(); u++)
    {
        INFO( "d "<<diag[u]);
        CHECK( fabs( y1[u] - 1./sqrt( diag[u])) < 1e-8);
    }
}

TEST_CASE( "InvSqrtODE factories")
{
    // The factories return by value; the result must survive the move into
    // make_directODESolve
    Container y0( diag.size(), 1.), y1( y0);
    unsigned number = 0;
    SECTION( "make_inv_sqrtodeCG")
    {
        auto sqrtinv = dg::mat::make_directODESolve(
                dg::mat::make_inv_sqrtodeCG( diagonal, 1., diag, 1e-12),
                std::string("Dormand-Prince-7-4-5"), 1e-10, 1e-12, number);
        sqrtinv( y0, y1);
    }
    SECTION( "make_inv_sqrtodeTri")
    {
        dg::TriDiagonal<Container> T( diag.size());
        for( unsigned u=0; u<diag.size(); u++)
        {
            T.M[u] = T.P[u] = 0.;
            T.O[u] = diag[u];
        }
        auto sqrtinv = dg::mat::make_directODESolve(
                dg::mat::make_inv_sqrtodeTri( T, diag),
                std::string("Dormand-Prince-7-4-5"), 1e-10, 1e-12, number);
        sqrtinv( y0, y1);
    }
    for( unsigned u=0; u<diag.size(); u++)
    {
        INFO( "d "<<diag[u]);
        CHECK( fabs( y1[u] - 1./sqrt( diag[u])) < 1e-8);
    }
}
