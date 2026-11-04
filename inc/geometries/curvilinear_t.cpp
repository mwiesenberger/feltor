#include <iostream>

#include "dg/algorithm.h"
#include "ribeiro.h"
#include "catch2/catch_all.hpp"

// Flux function sign*(1-r^2) with r^2 = (R-10)^2+Z^2
static dg::geo::CylindricalFunctorsLvl2 circular_psi( double sign)
{
    return dg::geo::CylindricalFunctorsLvl2(
        [=]( double R, double Z){ return sign*(1. - (R-10.)*(R-10.) - Z*Z);},
        [=]( double R, double  ){ return -2.*sign*(R-10.);},
        [=]( double  , double Z){ return -2.*sign*Z;},
        [=]( double  , double  ){ return -2.*sign;},
        [=]( double  , double  ){ return 0.;},
        [=]( double  , double  ){ return -2.*sign;});
}

// The metric of a generator must be the inverse Jacobian of its mapped points
static void check_metric( const dg::geo::aGenerator2d& generator, double zeta, double eta)
{
    const double hz = 1e-3*generator.width(), he = 1e-3;
    thrust::host_vector<double> zeta1d{ zeta-hz, zeta, zeta+hz}, eta1d{ eta-he, eta, eta+he};
    thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
    generator.generate( zeta1d, eta1d, x, y, zetaX, zetaY, etaX, etaY);
    const double xZeta = (x[5]-x[3])/2./hz, xEta = (x[7]-x[1])/2./he;
    const double yZeta = (y[5]-y[3])/2./hz, yEta = (y[7]-y[1])/2./he;
    const double det = xZeta*yEta - xEta*yZeta;
    const double eps = 1e-5;
    CHECK_THAT( zetaX[4], Catch::Matchers::WithinRel( +yEta/det, eps) ||
            Catch::Matchers::WithinAbs( +yEta/det, eps));
    CHECK_THAT( zetaY[4], Catch::Matchers::WithinRel( -xEta/det, eps) ||
            Catch::Matchers::WithinAbs( -xEta/det, eps));
    CHECK_THAT( etaX[4],  Catch::Matchers::WithinRel( -yZeta/det, eps) ||
            Catch::Matchers::WithinAbs( -yZeta/det, eps));
    CHECK_THAT( etaY[4],  Catch::Matchers::WithinRel( +xZeta/det, eps) ||
            Catch::Matchers::WithinAbs( +xZeta/det, eps));
}

// The flux function on the boundaries of the computational domain
static void check_boundaries( const dg::geo::aGenerator2d& generator,
    const dg::geo::CylindricalFunctorsLvl2& psi, double psi_0, double psi_1)
{
    const double w = generator.width();
    thrust::host_vector<double> zeta1d{ 1e-6*w, (1.-1e-6)*w}, eta1d{ 1.1};
    thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
    generator.generate( zeta1d, eta1d, x, y, zetaX, zetaY, etaX, etaY);
    const double eps = 1e-4;
    CHECK_THAT( psi.f()( x[0], y[0]), Catch::Matchers::WithinAbs( psi_0, eps));
    CHECK_THAT( psi.f()( x[1], y[1]), Catch::Matchers::WithinAbs( psi_1, eps));
}

// All parameters of the flux surfaces r^2 = 0.25 and r^2 = 0.75
TEST_CASE( "Ribeiro")
{
    // psi either falls or grows with r
    const double sign = GENERATE( +1., -1.);
    // from the inner to the outer flux surface or the other way around
    const bool outward = GENERATE( true, false);
    const int mode = GENERATE( 0, 1);
    INFO( "sign "<<sign<<" outward "<<outward<<" mode "<<mode);
    const auto psi = circular_psi( sign);
    const double psi_0 = sign*(outward ? 0.75 : 0.25);
    const double psi_1 = sign*(outward ? 0.25 : 0.75);
    const dg::geo::Ribeiro generator( psi, psi_0, psi_1, 10.1, 0., mode);
    // zeta = int |f| dpsi with |f| = 1/2/r^2 for mode 0 and |f| = 1/r for mode 1
    CHECK_THAT( generator.width(), Catch::Matchers::WithinRel(
        mode == 0 ? log(3.)/2. : sqrt(3.)-1., 1e-8));
    CHECK_THAT( generator.height(), Catch::Matchers::WithinRel( 2.*M_PI, 1e-14));
    check_metric( generator, 0.3*generator.width(), 1.1);
    check_boundaries( generator, psi, psi_0, psi_1);
}
