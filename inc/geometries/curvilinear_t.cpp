#include <iostream>

#include "dg/algorithm.h"
#include "curvilinearX.h"
#include "flux.h"
#include "hector.h"
#include "polar.h"
#include "refined_curvilinearX.h"
#include "ribeiro.h"
#include "separatrix_orthogonal.h"
#include "simple_orthogonal.h"
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

// The coordinate zeta is linear in psi
static void check_linear( const dg::geo::aGenerator2d& generator,
    const dg::geo::CylindricalFunctorsLvl2& psi, double psi_0, double psi_1)
{
    const double w = generator.width();
    thrust::host_vector<double> zeta1d{ 0.3*w}, eta1d{ 1.1};
    thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
    generator.generate( zeta1d, eta1d, x, y, zetaX, zetaY, etaX, etaY);
    CHECK_THAT( psi.f()( x[0], y[0]), Catch::Matchers::WithinAbs(
        psi_0 + 0.3*(psi_1-psi_0), 1e-10));
}

// grad zeta * grad eta = 0
static void check_orthogonal( const dg::geo::aGenerator2d& generator)
{
    CHECK( generator.isOrthogonal());
    thrust::host_vector<double> zeta1d{ 0.3*generator.width()}, eta1d{ 1.1};
    thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
    generator.generate( zeta1d, eta1d, x, y, zetaX, zetaY, etaX, etaY);
    CHECK_THAT( zetaX[0]*etaX[0] + zetaY[0]*etaY[0], Catch::Matchers::WithinAbs( 0, 1e-10));
}

// The eta coordinate closes after 2 Pi
static void check_period( const dg::geo::aGenerator2d& generator)
{
    const unsigned N = 64;
    thrust::host_vector<double> zeta1d{ 0.3*generator.width()}, eta1d(N+1);
    for( unsigned i=0; i<=N; i++)
        eta1d[i] = 2.*M_PI*i/N;
    thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
    generator.generate( zeta1d, eta1d, x, y, zetaX, zetaY, etaX, etaY);
    CHECK_THAT( x[N], Catch::Matchers::WithinAbs( x[0], 1e-10));
    CHECK_THAT( y[N], Catch::Matchers::WithinAbs( y[0], 1e-10));
}

TEST_CASE( "Polar generators")
{
    SECTION( "Polar")
    {
        const dg::geo::PolarGenerator generator( 1, 2);
        CHECK( generator.width() == 1);
        CHECK( generator.height() == 2.*M_PI);
        check_metric( generator, 0.3, 1.1);
        check_orthogonal( generator);
        check_period( generator);
    }
    SECTION( "LogPolar")
    {
        const dg::geo::LogPolarGenerator generator( 1, 2);
        CHECK_THAT( generator.width(), Catch::Matchers::WithinRel( log(2.), 1e-14));
        check_metric( generator, 0.3, 1.1);
        check_orthogonal( generator);
        check_period( generator);
    }
}

// Between the flux surfaces r^2 = 0.7 and r^2 = 0.75: a narrow domain keeps the
// expensive generators fast
TEST_CASE( "Generators on circular flux surfaces")
{
    // psi either falls or grows with r
    const double sign = GENERATE( +1., -1.);
    // from the inner to the outer flux surface or the other way around
    const bool outward = GENERATE( true, false);
    INFO( "sign "<<sign<<" outward "<<outward);
    const auto psi = circular_psi( sign);
    const double r0 = sqrt( outward ? 0.7 : 0.75);
    const double psi_0 = sign*(outward ? 0.3 : 0.25);
    const double psi_1 = sign*(outward ? 0.25 : 0.3);
    const double delta_psi = 0.05;
    // f_0 = 1/2/r^2 for mode 0 and 1/r for mode 1 on the first surface
    auto first_f0 = [=]( int mode){ return mode == 0 ? 0.5/r0/r0 : 1./r0;};
    SECTION( "Ribeiro")
    {
        const int mode = GENERATE( 0, 1);
        INFO( "mode "<<mode);
        const dg::geo::Ribeiro generator( psi, psi_0, psi_1, 10.1, 0., mode);
        // zeta = int |f| dpsi with |f| = 1/2/r^2 for mode 0 and |f| = 1/r for mode 1
        CHECK_THAT( generator.width(), Catch::Matchers::WithinRel(
            mode == 0 ? log(0.75/0.7)/2. : 2.*(sqrt(0.75)-sqrt(0.7)), 1e-8));
        CHECK_THAT( generator.height(), Catch::Matchers::WithinRel( 2.*M_PI, 1e-14));
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_boundaries( generator, psi, psi_0, psi_1);
        check_period( generator);
    }
    SECTION( "Hector")
    {
        const dg::geo::Hector<> generator( psi, psi_0, psi_1, 10.1, 0., 7, 2, 10, 1e-8);
        const dg::geo::aGenerator2d& base = generator;
        CHECK( generator.isConformal());
        // The conformal map of a ring has the width log(r_1/r_0)
        CHECK_THAT( base.width(), Catch::Matchers::WithinRel( log(0.75/0.7)/2., 1e-6));
        check_metric( base, 0.3*base.width(), 1.1);
        check_boundaries( base, psi, psi_0, psi_1);
        check_orthogonal( base);
        // Cauchy-Riemann
        thrust::host_vector<double> zeta1d{ 0.3*base.width()}, eta1d{ 1.1};
        thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
        base.generate( zeta1d, eta1d, x, y, zetaX, zetaY, etaX, etaY);
        CHECK_THAT( zetaX[0]*zetaX[0], Catch::Matchers::WithinRel( etaY[0]*etaY[0], 1e-6));
        CHECK_THAT( zetaY[0]*zetaY[0], Catch::Matchers::WithinRel( etaX[0]*etaX[0], 1e-6));
    }
    SECTION( "RibeiroFluxGenerator")
    {
        const int mode = GENERATE( 0, 1);
        INFO( "mode "<<mode);
        const double f0 = first_f0( mode);
        const dg::geo::RibeiroFluxGenerator generator( psi, psi_0, psi_1, 10.1, 0., mode);
        CHECK_THAT( generator.width(), Catch::Matchers::WithinRel( delta_psi*f0, 1e-8));
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_linear( generator, psi, psi_0, psi_1);
        check_period( generator);
    }
    SECTION( "SimpleOrthogonal")
    {
        const int mode = GENERATE( 0, 1);
        INFO( "mode "<<mode);
        const double f0 = first_f0( mode);
        const dg::geo::SimpleOrthogonal generator( psi, psi_0, psi_1, 10.1, 0., psi_0, mode);
        CHECK_THAT( generator.width(), Catch::Matchers::WithinRel( delta_psi*f0, 1e-8));
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_linear( generator, psi, psi_0, psi_1);
        check_orthogonal( generator);
        check_period( generator);
    }
    SECTION( "SimpleOrthogonal with first line between the boundaries")
    {
        const int mode = GENERATE( 0, 1);
        INFO( "mode "<<mode);
        const dg::geo::SimpleOrthogonal generator( psi, psi_0, psi_1, 10.1, 0.,
            0.5*(psi_0+psi_1), mode);
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_linear( generator, psi, psi_0, psi_1);
        check_orthogonal( generator);
    }
    SECTION( "SimpleOrthogonal with monitor metric")
    {
        const int mode = GENERATE( 0, 1);
        INFO( "mode "<<mode);
        const auto chi = dg::geo::make_LiseikinCollective( psi, 2., 0.5);
        const dg::geo::SimpleOrthogonal generator( psi, chi, psi_0, psi_1, 10.1, 0., psi_0, mode);
        CHECK( not generator.isOrthogonal());
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_linear( generator, psi, psi_0, psi_1);
    }
    SECTION( "FluxGenerator")
    {
        const int mode = GENERATE( 0, 1);
        INFO( "mode "<<mode);
        const double f0 = first_f0( mode);
        // f_0 = 1/q on the first surface for mode 0
        const double I0 = 20.;
        const double q = I0/2./sqrt( 100. - r0*r0);
        const dg::geo::CylindricalFunctorsLvl1 ipol( dg::geo::Constant(I0),
            dg::geo::Constant(0), dg::geo::Constant(0));
        const dg::geo::FluxGenerator generator( psi, ipol, psi_0, psi_1, 10.1, 0., mode);
        CHECK_THAT( generator.width(), Catch::Matchers::WithinRel(
            delta_psi*(mode == 0 ? 1./q : f0), 1e-8));
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_linear( generator, psi, psi_0, psi_1);
        check_period( generator);
        if( mode == 0)
        {
            // The Jacobian of the map is q R/I/f_0
            const double w = generator.width(), f = w/delta_psi;
            thrust::host_vector<double> zeta1d{ 0.3*w}, eta1d{ 1.1};
            thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
            generator.generate( zeta1d, eta1d, x, y, zetaX, zetaY, etaX, etaY);
            const double r2 = (x[0]-10.)*(x[0]-10.) + y[0]*y[0];
            const double qr = I0/2./sqrt( 100. - r2);
            CHECK_THAT( fabs( 1./(zetaX[0]*etaY[0] - zetaY[0]*etaX[0])),
                Catch::Matchers::WithinRel( qr*x[0]/I0/f, 1e-8));
        }
    }
}

TEST_CASE( "Hector with adaption")
{
    const auto psi = circular_psi( -1.);
    const double psi_0 = -0.75, psi_1 = -0.25;
    SECTION( "Orthogonal with weight function")
    {
        const auto chi = dg::geo::make_NablaPsiInvCollective( psi);
        const dg::geo::Hector<> generator( psi, chi, psi_0, psi_1, 10.1, 0., 7, 2, 10, 1e-6);
        CHECK( not generator.isConformal());
        check_orthogonal( generator);
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_boundaries( generator, psi, psi_0, psi_1);
    }
    SECTION( "Monitor metric")
    {
        const auto chi = dg::geo::make_LiseikinCollective( psi, 2., 0.5);
        const dg::geo::Hector<> generator( psi, chi, psi_0, psi_1, 10.1, 0., 7, 2, 10, 1e-6);
        CHECK( not generator.isConformal());
        CHECK( not generator.isOrthogonal());
        check_metric( generator, 0.3*generator.width(), 1.1);
        check_boundaries( generator, psi, psi_0, psi_1);
    }
}

TEST_CASE( "Curvilinear grids")
{
    SECTION( "Area of an annulus")
    {
        const dg::geo::PolarGenerator generator( 1, 2);
        const double area = 3.*M_PI;
        dg::geo::CurvilinearGrid2d grid( generator, 3, 8, 16, dg::DIR, dg::PER);
        const dg::HVec vol = dg::create::volume( grid);
        CHECK_THAT( dg::blas1::dot( vol, 1.), Catch::Matchers::WithinRel( area, 1e-10));
        // Also after changing the resolution
        grid.set( {4,4}, {6,12});
        const dg::HVec vol_set = dg::create::volume( grid);
        CHECK_THAT( dg::blas1::dot( vol_set, 1.), Catch::Matchers::WithinRel( area, 1e-10));
        CHECK( grid.bcx() == dg::DIR);
        CHECK( grid.bcy() == dg::PER);
    }
    SECTION( "Volume of a torus shell")
    {
        // Int R dR dZ dphi = 2 Pi R_0 Pi (r_1^2 - r_0^2)
        const auto psi = circular_psi( -1.);
        const dg::geo::SimpleOrthogonal generator( psi, -0.75, -0.25, 10.1, 0., -0.75, 1);
        const double exact = 2.*M_PI*10.*M_PI*0.5;
        dg::geo::CurvilinearProductGrid3d grid( generator, 3, 8, 32, 4);
        const dg::HVec vol = dg::create::volume( grid);
        CHECK_THAT( dg::blas1::dot( vol, 1.), Catch::Matchers::WithinRel( exact, 1e-8));
        const auto perp = grid.perp_grid();
        const dg::HVec area = dg::create::volume( *perp);
        CHECK_THAT( dg::blas1::dot( area, 1.), Catch::Matchers::WithinRel( M_PI*0.5, 1e-8));
        // Changing the number of planes does not change the volume
        grid.set( {3,3,1}, {8,32,8});
        const dg::HVec vol_set = dg::create::volume( grid);
        CHECK( vol_set.size() == 2*vol.size());
        CHECK_THAT( dg::blas1::dot( vol_set, 1.), Catch::Matchers::WithinRel( exact, 1e-8));
        // The metric is diagonal for orthogonal grids (value 0 is the zero vector)
        const auto metric = grid.metric();
        CHECK( metric.idx(0,1) == 0);
        CHECK( dg::blas1::dot( metric.value(0,1), metric.value(0,1)) == 0);
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

TEST_CASE( "X-point monitor metrics")
{
    dg::file::WrappedJsonValue js( dg::file::string2Json( geometry_params_Xpoint),
            dg::file::error::is_silent);
    const dg::geo::solovev::Parameters gp( js);
    const dg::geo::TokamakMagneticField mag = dg::geo::createSolovevField( gp);
    const auto psi = mag.get_psip();
    double R_X = gp.R_0-1.1*gp.triangularity*gp.a, Z_X = -1.1*gp.elongation*gp.a;
    auto check = [&]( const dg::geo::CylindricalSymmTensorLvl1& chi, double R, double Z){
        const double det = chi.xx()(R,Z)*chi.yy()(R,Z) - chi.xy()(R,Z)*chi.xy()(R,Z);
        CHECK_THAT( det, Catch::Matchers::WithinRel( 1., 1e-12));
        // The Laplacian div( chi grad psi) vanishes at the X-point
        const double lap = dg::geo::detail::LaplaceChiPsi( psi, chi)(R,Z);
        CHECK_THAT( lap, Catch::Matchers::WithinAbs( 0., 1e-14));
    };
    SECTION( "Constant")
    {
        const auto chi = dg::geo::make_Xconst_monitor( psi, R_X, Z_X);
        check( chi, R_X, Z_X);
    }
    SECTION( "Bump")
    {
        const auto chi = dg::geo::make_Xbump_monitor( psi, R_X, Z_X, 10., 10.);
        check( chi, R_X, Z_X);
        // Away from the X-point
        const double R = R_X+3, Z = Z_X+2, h = 1e-3;
        auto dx = [=]( const dg::geo::CylindricalFunctor& f){
            return (f(R+h,Z) - f(R-h,Z))/2./h;};
        auto dy = [=]( const dg::geo::CylindricalFunctor& f){
            return (f(R,Z+h) - f(R,Z-h))/2./h;};
        CHECK_THAT( chi.divX()(R,Z), Catch::Matchers::WithinRel( dx( chi.xx()) + dy( chi.xy()), 1e-5));
        CHECK_THAT( chi.divY()(R,Z), Catch::Matchers::WithinRel( dx( chi.xy()) + dy( chi.yy()), 1e-5));
    }
}

TEST_CASE( "SeparatrixOrthogonal")
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
    const double psi_0 = 0.3*psi.f()( R_O, Z_O), fx = 0.25, fy = 1./22.;
    const int mode = GENERATE( 0, 1);
    INFO( "mode "<<mode);
    const dg::geo::SeparatrixOrthogonal generator( psi, chi, psi_0, R_X, Z_X, R_O, Z_O, mode);
    SECTION( "Metric and flux function inside and outside of the separatrix")
    {
        const double zeta0 = generator.zeta0( fx), zeta1 = generator.zeta1( fx);
        CHECK_THAT( zeta1, Catch::Matchers::WithinRel( -fx/(1.-fx)*zeta0, 1e-14));
        for( double zeta : {0.5*zeta0, 0.5*zeta1})
        {
            INFO( "zeta "<<zeta);
            const double hz = 1e-3*fabs(zeta0), eta = 1.1, he = 1e-3;
            thrust::host_vector<double> zeta1d{ zeta-hz, zeta, zeta+hz}, eta1d{ eta-he, eta, eta+he};
            thrust::host_vector<double> x, y, zetaX, zetaY, etaX, etaY;
            generator.generate( zeta1d, eta1d, 0, 3, x, y, zetaX, zetaY, etaX, etaY);
            const double xZeta = (x[5]-x[3])/2./hz, xEta = (x[7]-x[1])/2./he;
            const double yZeta = (y[5]-y[3])/2./hz, yEta = (y[7]-y[1])/2./he;
            const double det = xZeta*yEta - xEta*yZeta;
            const double eps = 1e-5;
            CHECK_THAT( zetaX[4], Catch::Matchers::WithinRel( +yEta/det, eps));
            CHECK_THAT( zetaY[4], Catch::Matchers::WithinRel( -xEta/det, eps));
            CHECK_THAT( etaX[4],  Catch::Matchers::WithinRel( -yZeta/det, eps));
            CHECK_THAT( etaY[4],  Catch::Matchers::WithinRel( +xZeta/det, eps));
            // zeta is linear in psi and zero on the separatrix
            const double eps_psi = mode == 0 ? 2e-3 : 1e-6;
            CHECK_THAT( psi.f()(x[4],y[4]), Catch::Matchers::WithinAbs( psi_0*zeta/zeta0, eps_psi));
        }
    }
    SECTION( "Refined grid covers the same area")
    {
        const dg::geo::CurvilinearGridX2d grid( generator, fx, fy, 3, 8, 44);
        const dg::HVec vol = dg::create::volume( grid);
        const dg::EquidistXRefinement refine( 0, 2, 1, 1);
        const dg::geo::CurvilinearRefinedGridX2d refined( refine, generator, fx, fy, 3, 8, 44);
        const dg::HVec volr = dg::create::volume( refined);
        CHECK_THAT( dg::blas1::dot( volr, 1.), Catch::Matchers::WithinRel(
            dg::blas1::dot( vol, 1.), 1e-2));
    }
}

TEST_CASE( "X-point grids with default boundary conditions")
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
    const dg::geo::SeparatrixOrthogonal generator( psi, chi, 0.3*psi.f()( R_O, Z_O),
            R_X, Z_X, R_O, Z_O, 1);
    const dg::geo::CurvilinearGridX2d grid2d( generator, 0.25, 1./22., 3, 4, 22);
    CHECK( grid2d.bcx() == dg::DIR);
    CHECK( grid2d.bcy() == dg::NEU);
    const dg::geo::CurvilinearProductGridX3d grid3d( generator, 0.25, 1./22., 3, 4, 22, 2);
    CHECK( grid3d.bcy() == dg::NEU);
    CHECK( grid3d.bcz() == dg::PER);
}
