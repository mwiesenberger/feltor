#include <iostream>

#include "dg/algorithm.h"
#include "average.h"
#include "toroidal.h"
#include "catch2/catch_all.hpp"

// Circular flux surfaces psi = 1 - r^2 with r^2 = (R-R_0)^2 + Z^2 and I = I_0
// In polar coordinates all integrals are analytic
TEST_CASE( "Flux surface integrals")
{
    const double R_0 = 10, I_0 = 20;
    const dg::geo::TokamakMagneticField mag = dg::geo::createCircularField( R_0, I_0, 1, 1);
    const dg::Grid2d grid( R_0-1.5, R_0+1.5, -1.5, 1.5, 3, 200, 200);
    const dg::HVec one = dg::evaluate( dg::one, grid);
    // 1 for Z > 0, cuts away half of each flux surface
    const dg::HVec upper = dg::evaluate( []( double, double Z){ return Z>0 ? 1. : 0.;}, grid);
    const double psi0 = GENERATE( 0.25, 0.5, 0.75);
    INFO( "psi0 "<<psi0);
    const double r2 = 1-psi0;

    SECTION( "FluxSurfaceIntegral")
    {
        // int r dr dtheta delta(1 - r^2 - psi0) = pi
        dg::geo::FluxSurfaceIntegral<dg::HVec> fsi( grid, mag);
        CHECK_THAT( fsi( psi0), Catch::Matchers::WithinRel( M_PI, 1e-5));
        fsi.set_right( upper);
        CHECK_THAT( fsi( psi0), Catch::Matchers::WithinRel( M_PI/2., 1e-3));
        // A wider numerical delta function
        dg::geo::FluxSurfaceIntegral<dg::HVec> wide( grid, mag, 2.);
        CHECK( wide.get_deltapsi() == 2.*fsi.get_deltapsi());
        CHECK_THAT( wide( psi0), Catch::Matchers::WithinRel( M_PI, 1e-5));
    }
    SECTION( "FluxVolumeIntegral")
    {
        // psi <= psi0 is the box without the disc of radius r
        dg::geo::FluxVolumeIntegral<dg::HVec> fvi( (dg::CartesianGrid2d)grid, mag);
        CHECK_THAT( fvi( psi0), Catch::Matchers::WithinRel( 9. - M_PI*r2, 1e-3));
        fvi.set_left( upper);
        CHECK_THAT( fvi( psi0), Catch::Matchers::WithinRel( 4.5 - M_PI*r2/2., 1e-3));
    }
    SECTION( "FluxSurfaceAverage")
    {
        // The average of R over a circle with uniform measure is R_0
        const dg::HVec R = dg::evaluate( dg::cooX2d, grid);
        dg::geo::FluxSurfaceAverage<dg::HVec> fsa( grid, mag, R, one);
        CHECK_THAT( fsa( psi0), Catch::Matchers::WithinRel( R_0, 1e-5));
        fsa.set_container( one);
        CHECK_THAT( fsa( psi0), Catch::Matchers::WithinRel( 1., 1e-6));
        dg::geo::FluxSurfaceIntegral<dg::HVec> fsi( grid, mag);
        CHECK( fsa.get_deltapsi() == fsi.get_deltapsi());
    }
    SECTION( "Safety factor")
    {
        // q = I_0/2/sqrt( R_0^2 - r^2)
        const double q = I_0/2./sqrt(R_0*R_0 - r2);
        dg::geo::SafetyFactor direct( mag);
        CHECK_THAT( direct( psi0), Catch::Matchers::WithinRel( q, 1e-6));
        dg::geo::SafetyFactorAverage average( grid, mag);
        CHECK_THAT( average( psi0), Catch::Matchers::WithinRel( q, 1e-5));
        // The sign follows the current
        dg::geo::SafetyFactor negative( dg::geo::createCircularField( R_0, -I_0, 1, 1));
        CHECK_THAT( negative( psi0), Catch::Matchers::WithinRel( -q, 1e-6));
        dg::geo::SafetyFactorAverage half( grid, mag);
        half.set_weights( upper);
        CHECK_THAT( half( psi0), Catch::Matchers::WithinRel( q/2., 1e-3));
    }
}
