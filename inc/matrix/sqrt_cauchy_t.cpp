#include <iostream>

#include "sqrt_cauchy.h"

#include "catch2/catch_all.hpp"

using Container = dg::HVec;

// Multiplication with a diagonal d has analytic matrix functions
static const std::vector<double> diag_values = {0.5, 1., 2., 3., 10.};
static const Container diag( diag_values.begin(), diag_values.end());
static void diagonal( const Container& y, Container& z)
{
    dg::blas1::pointwiseDot( diag, y, z);
}

static dg::mat::DirectSqrtCauchy<Container> make_direct( int exp)
{
    return { diagonal, diag, 1e-12, 40, {0.5, 10.}, exp};
}

TEST_CASE( "DirectSqrtCauchy")
{
    // x = d^{exp/2} b
    int exp = GENERATE( +1, -1);
    INFO( "exp "<<exp);
    Container b( diag.size(), 1.), x( b);
    dg::mat::DirectSqrtCauchy<Container> direct;
    SECTION( "Constructor")
    {
        direct = make_direct( exp);
    }
    SECTION( "construct")
    {
        direct.construct( diagonal, diag, 1e-12, 40, std::array<double,2>{0.5, 10.}, exp);
    }
    SECTION( "Copy")
    {
        auto other = make_direct( exp);
        direct = other;
    }
    direct( b, x);
    for( unsigned u=0; u<diag.size(); u++)
    {
        INFO( "d "<<diag[u]);
        CHECK( fabs( x[u] - pow( diag[u], exp/2.)) < 1e-10);
    }
}
