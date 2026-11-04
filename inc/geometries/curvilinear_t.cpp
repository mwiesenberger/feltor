#include <iostream>

#include "dg/algorithm.h"
#include "flux.h"
#include "hector.h"
#include "ribeiro.h"
#include "separatrix_orthogonal.h"
#include "solovev.h"
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

TEST_CASE( "Hector")
{
    const double sign = GENERATE( +1., -1.);
    const bool outward = GENERATE( true, false);
    INFO( "sign "<<sign<<" outward "<<outward);
    const auto psi = circular_psi( sign);
    const double psi_0 = sign*(outward ? 0.75 : 0.25);
    const double psi_1 = sign*(outward ? 0.25 : 0.75);
    const dg::geo::Hector<> generator( psi, psi_0, psi_1, 10.1, 0., 7, 2, 10, 1e-8);
    SECTION( "The conformal map of a ring has the width log(r_1/r_0)")
    {
        const dg::geo::aGenerator2d& base = generator;
        CHECK( generator.isConformal());
        CHECK( base.isOrthogonal());
        CHECK_THAT( base.width(), Catch::Matchers::WithinRel( log(3.)/2., 1e-6));
        check_metric( base, 0.3*base.width(), 1.1);
        check_boundaries( base, psi, psi_0, psi_1);
        thrust::host_vector<double> zeta1d{ 0.3*base.width()}, eta1d{ 1.1};
        thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
        base.generate( zeta1d, eta1d, x, y, zetaX, zetaY, etaX, etaY);
        // Cauchy-Riemann
        CHECK_THAT( zetaX[0]*zetaX[0], Catch::Matchers::WithinRel( etaY[0]*etaY[0], 1e-6));
        CHECK_THAT( zetaY[0]*zetaY[0], Catch::Matchers::WithinRel( etaX[0]*etaX[0], 1e-6));
    }
}

TEST_CASE( "Flux surface with psi = 0")
{
    // The surface r = 1 has psi = 0 and lies in the middle of the domain
    const double psi_0 = 0.5, psi_1 = -0.5;
    const auto psi = circular_psi( 1.);
    const dg::geo::CylindricalFunctorsLvl1 ipol( dg::geo::Constant(20),
            dg::geo::Constant(0), dg::geo::Constant(0));
    SECTION( "FluxGenerator")
    {
        const int mode = GENERATE( 0, 1);
        INFO( "mode "<<mode);
        const dg::geo::FluxGenerator generator( psi, ipol, psi_0, psi_1, 10.1, 0., mode);
        check_metric( generator, 0.5*generator.width(), 1.1);
    }
    SECTION( "RibeiroFluxGenerator")
    {
        const int mode = GENERATE( 0, 1);
        INFO( "mode "<<mode);
        const dg::geo::RibeiroFluxGenerator generator( psi, psi_0, psi_1, 10.1, 0., mode);
        check_metric( generator, 0.5*generator.width(), 1.1);
    }
}

const std::string geometry_params_Xpoint = R"asdf({
    "A" : 0.0,
    "PP": 1,
    "PI": 1,
    "c" :[  0.07350114445500399706283007092406934834526,
           -0.08662417436317227513877947632069712210813,
           -0.1463931543401102620740934776490506239925,
           -0.07631237100536276213126232216649739043965,
            0.09031790113794227394476271394334515457567,
           -0.09157541239018724584036670247895160625891,
           -0.003892282979837564486424586266476650443202,
            0.04271891225076417603805495295590637082745,
            0.2275545646002791311716859154040182853650,
           -0.1304724136017769544849838714185080700328,
           -0.03006974108476955225335835678782676287818,
            0.004212671892103931178531621087067962015783 ],
    "R_0"                : 547.891714877869,
    "inverseaspectratio" : 0.41071428571428575,
    "elongation"         : 1.75,
    "triangularity"      : 0.47,
    "equilibrium"  : "solovev",
    "description" : "standardX"
})asdf";

TEST_CASE( "SeparatrixOrthogonalAdaptor")
{
    dg::file::WrappedJsonValue js( dg::file::string2Json( geometry_params_Xpoint),
            dg::file::error::is_silent);
    const dg::geo::solovev::Parameters gp( js);
    const dg::geo::TokamakMagneticField mag = dg::geo::createSolovevField( gp);
    const auto psi = mag.get_psip();
    double R_X = gp.R_0-1.1*gp.triangularity*gp.a, Z_X = -1.1*gp.elongation*gp.a;
    double R_O = gp.R_0, Z_O = 0.;
    dg::geo::findOpoint( psi, R_O, Z_O);
    const auto chi = dg::geo::make_Xconst_monitor( psi, R_X, Z_X);
    const double psi_0 = 0.3*psi.f()( R_O, Z_O);
    const unsigned mode = 1;
    SECTION( "The default ends on the separatrix")
    {
        const dg::geo::SeparatrixOrthogonalAdaptor generator( psi, chi, psi_0, R_X,
                Z_X, R_O, Z_O, mode);
        CHECK_THAT( generator.height(), Catch::Matchers::WithinRel( 2.*M_PI, 1e-14));
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_boundaries( generator, psi, psi_0, 0.);
    }
    SECTION( "A fraction of the points in the scrape-off layer")
    {
        const double fx = 0.25;
        const dg::geo::SeparatrixOrthogonalAdaptor generator( psi, chi, psi_0, R_X,
                Z_X, R_O, Z_O, mode, false, fx);
        const double w = generator.width();
        thrust::host_vector<double> eta1d{ 1.1}, x, y, zetaX, zetaY, etaX, etaY;
        // The separatrix is at 3/4 w, psi is linear in zeta
        generator.generate( {0.25*w, 0.5*w, 0.75*w, w}, eta1d, x, y, zetaX, zetaY, etaX, etaY);
        const double psi_1 = -fx/(1.-fx)*psi_0;
        CHECK_THAT( psi.f()( x[2], y[2]), Catch::Matchers::WithinAbs( 0., 1e-4));
        CHECK_THAT( psi.f()( x[3], y[3]), Catch::Matchers::WithinAbs( psi_1, 1e-4));
        thrust::host_vector<double> three{ 0.25*w, 0.5*w, 0.75*w};
        CHECK_THROWS_AS( generator.generate( three, eta1d, x, y, zetaX, zetaY, etaX, etaY),
                dg::Error);
    }
}
