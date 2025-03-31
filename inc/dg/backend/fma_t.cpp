#include <iostream>
#include <complex>
#include "config.h"

#include "catch2/catch_all.hpp"

TEST_CASE( "Test FMA support")
{
    double t = 3, f = 4, ff = 5;
    CHECK( DG_FMA( t,f,ff) == 17);
    SECTION( "std")
    {
    std::complex<double> z0( 1,2), z1( 3,4);
    CHECK( DG_FMA( t, z0, z1) == std::complex<double>( 6, 10));
    std::complex<double> z2( 5,6), z3;
    z3 = z0*z1+z2;
    CHECK( DG_FMA( z0, z1, z2) == z3);
    z3 = 2*t+z2;
    CHECK( DG_FMA( 2.,  t, z2) == z3);
    z3 = t*z1+z2;
    CHECK( DG_FMA(  t, z1, z2) == z3);
    z3 = z1*t+z2;
    CHECK( DG_FMA(  z1, t, z2) == z3);
    }
    SECTION( "thrust")
    {
    thrust::complex<double> z0( 1,2), z1( 3,4);
    CHECK( DG_FMA( t, z0, z1) == thrust::complex<double>( 6, 10));
    thrust::complex<double> z2( 5,6), z3;
    z3 = z0*z1+z2;
    CHECK( DG_FMA( z0, z1, z2) == z3);
    z3 = 2*t+z2;
    CHECK( DG_FMA( 2.,  t, z2) == z3);
    z3 = t*z1+z2;
    CHECK( DG_FMA(  t, z1, z2) == z3);
    z3 = z1*t+z2;
    CHECK( DG_FMA(  z1, t, z2) == z3);
    }
}
