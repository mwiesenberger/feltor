#include <iostream>

#include "dg/algorithm.h"
#include "dg/file/json_utilities.h"
#include "make_field.h"
#include "sheath.h"
#include "catch2/catch_all.hpp"

// In the field v = (0, 0.5, sigma) in (R,Z,phi) field lines are straight:
// dZ/dphi = 0.5/sigma and ds/dphi = 1/sigma. Along the field Z always increases
TEST_CASE( "Wall fieldline functors")
{
    const double sigma = GENERATE( +1., -1.);
    INFO( "sigma "<<sigma);
    const dg::geo::CylindricalVectorLvl0 vec( dg::geo::Constant(0),
            dg::geo::Constant(0.5), dg::geo::Constant(sigma));
    const dg::Grid2d box( 2, 4, -1, 1, 1, 1, 1);
    const double R = 3, eps = 1e-5;
    const std::string type = GENERATE( "phi", "s");
    INFO( "type "<<type);

    SECTION( "Distance to wall in both directions")
    {
        dg::geo::WallFieldlineDistance plus( vec, box, +10, 1e-6, type);
        dg::geo::WallFieldlineDistance minus( vec, box, -10, 1e-6, type);
        CHECK_THAT( plus( R, 0.2), Catch::Matchers::WithinAbs( +1.6, eps));
        CHECK_THAT( minus( R, 0.2), Catch::Matchers::WithinAbs( -2.4, eps));
    }
    SECTION( "Distance if wall is not reached")
    {
        dg::geo::WallFieldlineDistance plus( vec, box, +1, 1e-6, type);
        dg::geo::WallFieldlineDistance minus( vec, box, -1, 1e-6, type);
        CHECK_THAT( plus( R, 0.2), Catch::Matchers::WithinAbs( +1, eps));
        CHECK_THAT( minus( R, 0.2), Catch::Matchers::WithinAbs( -1, eps));
    }
    SECTION( "Distance if predicate is false")
    {
        const double factor = type == "phi" ? 1. : 1e6;
        dg::geo::WallFieldlineDistance plus( vec, box, +10, 1e-6, type,
                dg::geo::mod::nowhere);
        dg::geo::WallFieldlineDistance minus( vec, box, -10, 1e-6, type,
                dg::geo::mod::nowhere);
        CHECK( plus( R, 0.2) == +10*factor);
        CHECK( minus( R, 0.2) == -10*factor);
    }
    SECTION( "Coordinate is independent of the sign of maxPhi")
    {
        const double maxPhi = GENERATE( +10., -10.);
        INFO( "maxPhi "<<maxPhi);
        dg::geo::WallFieldlineCoordinate coord( vec, box, maxPhi, 1e-6, type);
        // In the box Z in [-1,1] the coordinate is simply Z
        for( double Z : {-0.9, -0.3, 0.2, 0.7})
        {
            INFO( "Z "<<Z);
            CHECK_THAT( coord( R, Z), Catch::Matchers::WithinAbs( Z, eps));
        }
    }
    SECTION( "Coordinate if the wall is not reached")
    {
        // No wall reached
        dg::geo::WallFieldlineCoordinate coord( vec, box, 1, 1e-6, type);
        CHECK_THAT( coord( R, 0.2), Catch::Matchers::WithinAbs( 0, eps));
        // Only the wall in field direction reached: (1-0.2)/1.2
        CHECK_THAT( coord( R, +0.9), Catch::Matchers::WithinAbs( +2./3., eps));
        // Only the wall against field direction reached
        CHECK_THAT( coord( R, -0.9), Catch::Matchers::WithinAbs( -2./3., eps));
    }
    SECTION( "Coordinate if predicate is false")
    {
        dg::geo::WallFieldlineCoordinate coord( vec, box, 10, 1e-6, type,
                dg::geo::mod::nowhere);
        CHECK( coord( R, 0.2) == 0);
    }
}

TEST_CASE( "Sheath region")
{
    // The closed field line region is r = sqrt( (R-3)^2 + Z^2) < 1
    const double I0 = GENERATE( +10., -10.);
    INFO( "I_0 "<<I0);
    dg::geo::TokamakMagneticField mag = dg::geo::createCircularField( 3, I0, 1, 1);
    dg::file::WrappedJsonValue jsmod( dg::file::string2Json(
            R"({"boundary":0.0625,"alpha":0.015625})"), dg::file::error::is_silent);
    const dg::Grid2d box( 1.5, 4.5, -1.5, 1.5, 1, 1, 1);
    dg::Grid2d walls( box);
    dg::geo::CylindricalFunctor sheath;
    SECTION( "No wall")
    {
        dg::geo::createSheathRegion( jsmod, mag, dg::geo::Constant(0), walls, sheath);
        CHECK( walls.x0() == box.x0());
        CHECK( walls.x1() == box.x1());
        CHECK( walls.y0() == box.y0());
        CHECK( walls.y1() == box.y1());
        CHECK( sheath( 3.2, 0.0) == 0); // closed field line
        // field line is a circle with r = 1.57 and reaches R=4.5 after |phi| = 0.11
        CHECK( sheath( 4.45, 0.6) == 1);
        // field line is a circle with r = 1.48 and never leaves the box
        CHECK( sheath( 4.45, 0.3) == 0);
        CHECK( sheath( 4.2, 0.0) == 0);
    }
    SECTION( "Boundaries in the wall are not sheath")
    {
        dg::geo::CylindricalFunctor wall = []( double R, double Z){
            return R > 4.4 || Z > 1.4 ? 1. : 0.;
        };
        dg::geo::createSheathRegion( jsmod, mag, wall, walls, sheath);
        CHECK( walls.x0() == box.x0());
        CHECK( walls.x1() == 1e10);
        CHECK( walls.y0() == box.y0());
        CHECK( walls.y1() == 1e10);
        CHECK( sheath( 4.45, 0.6) == 0); // in the wall
    }
}
