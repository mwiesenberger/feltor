
#include <iostream>
#include <iomanip>

#include "dg/algorithm.h"

#include "solovev.h"
#include "magnetic_field.h"
#include "catch2/catch_all.hpp"

// We correct PI so that B(R_O, 0) = 1 exactly for A = 0
const std::string geometry_params_Xpoint = R"asdf({
    "A": 0,
    "R_0": 0.545,
    "PP": 1,
    "PI": 1.042138271304694,
    "c":
    [
        0.085858385692821789,
        0.25157608973291837,
        -0.38061591778794789,
        -0.10218614802095409,
        0.36709261273307433,
        -0.39303943992277839,
        -0.016931807751497822,
        0,
        0,
        0,
        0,
        0
    ],
    "equilibrium": "solovev",
    "description" : "doubleX",
    "inverseaspectratio": 0.3211009174311926,
    "triangularity": 0.3,
    "elongation": 1.44
})asdf";

TEST_CASE("Magnetic field")
{
    dg::file::WrappedJsonValue js( dg::file::error::is_warning);
    js.asJson() = dg::file::string2Json( geometry_params_Xpoint);
    // A != 0 makes the current I depend on psi
    double A = GENERATE( 0., 1.);
    js.asJson()["A"] = A;

    dg::geo::solovev::Parameters gp(js);
    int sign = GENERATE( +1, -1);
    gp.pi = sign*gp.pi;
    INFO( "A "<<A<<" sign "<<sign);

    dg::geo::TokamakMagneticField mag = dg::geo::createSolovevField(gp);

    double R_O = gp.R_0;
    double Z_O = 0.;
    dg::geo::findOpoint( mag.get_psip(), R_O, Z_O);
    // The poloidal field vanishes at the O-point
    const double B_O = gp.R_0/R_O*fabs( mag.ipol()( R_O, Z_O));
    SECTION( "O-point")
    {
        CHECK( fabs(R_O - ( A == 0 ? 0.567965 : 0.569601)) < 1e-6);
        CHECK( fabs(Z_O ) < 1e-6);
        CHECK( fabs( mag.psipR()( R_O, Z_O)) < 1e-12);
        CHECK( fabs( mag.psipZ()( R_O, Z_O)) < 1e-12);
        if( A == 0)
            CHECK( fabs( B_O - 1) < 1e-12);
    }

    SECTION( "Toroidal field")
    {
        // Magnetic field strength at O-point is purely toroidal
        dg::geo::Bmodule bmod( mag);
        dg::geo::Btor btor(mag);
        dg::geo::InvB invB(mag);
        dg::geo::InvBtor invBtor(mag);
        dg::geo::LnB lnB(mag);
        CHECK ( fabs( bmod( R_O, 0) - B_O ) < 1e-12);
        CHECK ( fabs( btor( R_O, 0) - sign*B_O ) < 1e-12);
        CHECK ( fabs( invB( R_O, 0) - 1./B_O ) < 1e-12);
        CHECK ( fabs( invBtor( R_O, 0) - sign/B_O ) < 1e-12);
        CHECK ( fabs( lnB(  R_O, 0) - log( B_O) ) < 1e-12);

        double R = R_O + 0.2, Z = -0.1;
        CHECK( fabs( 1./bmod(R,Z) - invB(R,Z)) < 1e-12);
        CHECK( fabs( 1./btor(R,Z) - invBtor(R,Z)) < 1e-12);

    }
    SECTION( "Curvature")
    {
        // Curvatures at O-point are the same for all
        dg::geo::CurvatureNablaBR  curvKBR( mag, sign);
        dg::geo::CurvatureNablaBZ  curvKBZ( mag, sign);
        dg::geo::CurvatureKappaR   curvKKR( mag, sign);
        dg::geo::CurvatureKappaZ   curvKKZ( mag, sign);
        dg::geo::DivCurvatureKappa divCurvKappa( mag, sign);

        dg::geo::ToroidalCurvatureNablaBR  torcurvKBR( mag);
        dg::geo::ToroidalCurvatureNablaBZ  torcurvKBZ( mag);
        dg::geo::ToroidalCurvatureKappaR   torcurvKKR( mag);
        dg::geo::ToroidalCurvatureKappaZ   torcurvKKZ( mag);
        dg::geo::ToroidalDivCurvatureKappa tordivCurvKappa( mag);

        dg::geo::TrueCurvatureNablaBR  truecurvKBR( mag);
        dg::geo::TrueCurvatureNablaBZ  truecurvKBZ( mag);
        dg::geo::TrueCurvatureKappaR   truecurvKKR( mag);
        dg::geo::TrueCurvatureKappaZ   truecurvKKZ( mag);
        dg::geo::TrueDivCurvatureKappa truedivCurvKappa( mag);

        CHECK( fabs(curvKBR( R_O, Z_O) - torcurvKBR( R_O, Z_O)) < 1e-12);
        CHECK( fabs(curvKBZ( R_O, Z_O) - torcurvKBZ( R_O, Z_O)) < 1e-12);
        CHECK( fabs(curvKKR( R_O, Z_O) - torcurvKKR( R_O, Z_O)) < 1e-12);
        CHECK( fabs(curvKKZ( R_O, Z_O) - torcurvKKZ( R_O, Z_O)) < 1e-12);
        CHECK( fabs(divCurvKappa( R_O, Z_O) - tordivCurvKappa( R_O, Z_O)) < 1e-12);

        CHECK( fabs(curvKBR( R_O, Z_O) - truecurvKBR( R_O, Z_O)) < 1e-12);
        CHECK( fabs(curvKBZ( R_O, Z_O) - truecurvKBZ( R_O, Z_O)) < 1e-12);
        CHECK( fabs(curvKKR( R_O, Z_O) - truecurvKKR( R_O, Z_O)) < 1e-12);
        CHECK( fabs(curvKKZ( R_O, Z_O) - truecurvKKZ( R_O, Z_O)) < 1e-12);
        CHECK( fabs(divCurvKappa( R_O, Z_O) - truedivCurvKappa( R_O, Z_O)) < 1e-12);

    }
    SECTION( "BHat")
    {
        auto bhat = dg::geo::createBHat( mag);
        auto torbhat = dg::geo::createToroidalBHat( mag);
        double R = R_O + 0.1, Z = 0.1;

        CHECK( fabs(bhat.x()(R,Z)/bhat.z()(R,Z)   - R*torbhat.x()( R,Z)) < 1e-12);
        CHECK( fabs(bhat.y()(R,Z)/bhat.z()(R,Z)   - R*torbhat.y()( R,Z)) < 1e-12);
        CHECK( fabs(1./R   - torbhat.z()( R,Z)) < 1e-12);
        CHECK ( fabs( torbhat.div()(R,Z) - ( bhat.divvvz()(R,Z)/R - bhat.x()(R,Z)/bhat.z()(R,Z)/R/R )) < 1e-12);
    }

    SECTION( "Divb")
    {
        dg::geo::Divb div(mag);
        dg::geo::GradLnB gradLnB(mag);
        dg::geo::ToroidalDivb tordiv(mag);
        dg::geo::ToroidalGradLnB torgradLnB(mag);
        dg::geo::DivVVP divvvp(mag);
        CHECK( fabs(div( R_O, 0)  - tordiv( R_O,0))< 1e-12);
        double R = R_O + 0.1, Z = 0.1;
        CHECK( fabs( div( R,Z) + gradLnB(  R,Z)) < 1e-12);

        CHECK( fabs( tordiv( R,Z) + torgradLnB(  R,Z)) < 1e-12);

        dg::geo::BFieldR BR(mag);
        dg::geo::BFieldP Bphi(mag);
        CHECK ( fabs( tordiv(R,Z) - ( divvvp(R,Z)/R - BR(R,Z)/Bphi(R,Z)/R/R )) < 1e-12);

        if( A != 0)
        {
            REQUIRE( mag.ipolR()(R,Z) != 0);
            REQUIRE( mag.ipolZ()(R,Z) != 0);
        }
        // 4th order central differences
        const double h = 1e-3;
        auto dR = [h]( const dg::geo::CylindricalFunctor& f, double r, double z){
            return (-f(r+2*h,z)+8*f(r+h,z)-8*f(r-h,z)+f(r-2*h,z))/(12*h);
        };
        auto dZ = [h]( const dg::geo::CylindricalFunctor& f, double r, double z){
            return (-f(r,z+2*h)+8*f(r,z+h)-8*f(r,z-h)+f(r,z-2*h))/(12*h);
        };
        dg::geo::CylindricalFunctor bR = dg::geo::ToroidalBHatR(mag);
        dg::geo::CylindricalFunctor bZ = dg::geo::ToroidalBHatZ(mag);
        dg::geo::CylindricalFunctor Rb = [&bR]( double r, double z){ return r*bR(r,z);};
        dg::geo::Btor btor(mag);
        dg::geo::CylindricalFunctor lnB = [&btor]( double r, double z){
            return log( fabs( btor(r,z)));
        };
        const double eps = 1e-6;
        CHECK_THAT( torgradLnB(R,Z), Catch::Matchers::WithinRel(
                bR(R,Z)*dR( lnB, R, Z) + bZ(R,Z)*dZ( lnB, R, Z), eps));
        CHECK_THAT( tordiv(R,Z), Catch::Matchers::WithinRel(
                dR( Rb, R, Z)/R + dZ( bZ, R, Z), eps));
    }

    SECTION( "BR BZ")
    {
        dg::geo::BR bR(mag);
        dg::geo::BZ bZ(mag);

        dg::geo::ToroidalBR torbR(mag);
        dg::geo::ToroidalBZ torbZ(mag);

        CHECK( fabs(bR( R_O, Z_O) - sign*torbR( R_O, Z_O)) < 1e-12);
        CHECK( fabs(bZ( R_O, Z_O) - sign*torbZ( R_O, Z_O)) < 1e-12);
    }
}

// f = a (R-1)^p + b (Z-2)^2 with an extremum or saddle point at (1,2)
static dg::geo::CylindricalFunctorsLvl2 create_critical( double a, int p, double b)
{
    return dg::geo::CylindricalFunctorsLvl2(
        [=]( double R, double Z){ return a*pow(R-1,p) + b*(Z-2)*(Z-2);},
        [=]( double R, double  ){ return a*p*pow(R-1,p-1);},
        [=]( double  , double Z){ return 2*b*(Z-2);},
        [=]( double R, double  ){ return a*p*(p-1)*pow(R-1,p-2);},
        []( double, double){ return 0.;},
        [=]( double, double){ return 2*b;});
}

TEST_CASE("Critical points")
{
    double R = 1.5, Z = 2.5;
    SECTION( "Local minimum")
    {
        CHECK( dg::geo::findCriticalPoint( create_critical( 1,2,1), R, Z) == 1);
        CHECK( fabs( R-1) < 1e-14);
        CHECK( fabs( Z-2) < 1e-14);
        CHECK( dg::geo::findOpoint( create_critical( 1,2,1), R, Z) == 1);
        CHECK_THROWS_AS( dg::geo::findXpoint( create_critical( 1,2,1), R, Z), dg::Error);
    }
    SECTION( "Local maximum")
    {
        CHECK( dg::geo::findCriticalPoint( create_critical( -1,2,-1), R, Z) == 2);
        CHECK( fabs( R-1) < 1e-14);
        CHECK( fabs( Z-2) < 1e-14);
        CHECK( dg::geo::findOpoint( create_critical( -1,2,-1), R, Z) == 2);
    }
    SECTION( "Saddle point")
    {
        CHECK( dg::geo::findCriticalPoint( create_critical( 1,2,-1), R, Z) == 3);
        CHECK( fabs( R-1) < 1e-14);
        CHECK( fabs( Z-2) < 1e-14);
        CHECK_NOTHROW( dg::geo::findXpoint( create_critical( 1,2,-1), R, Z));
        CHECK_THROWS_AS( dg::geo::findOpoint( create_critical( 1,2,-1), R, Z), dg::Error);
    }
    SECTION( "No critical point")
    {
        // Hessian vanishes everywhere; input stays unchanged
        dg::geo::CylindricalFunctorsLvl2 constant( dg::geo::Constant(1),
                dg::geo::Constant(0), dg::geo::Constant(0), dg::geo::Constant(0),
                dg::geo::Constant(0), dg::geo::Constant(0));
        CHECK( dg::geo::findCriticalPoint( constant, R, Z) == 0);
        CHECK( R == 1.5);
        CHECK( Z == 2.5);
    }
    SECTION( "Convergence at a degenerate critical point")
    {
        // The Hessian vanishes at the critical point so Newton only converges
        // linearly with factor 2/3. The error is then 2 times the last step,
        // which must be below 1e-10 and not only below 1e-7
        CHECK( dg::geo::findCriticalPoint( create_critical( 1,4,1), R, Z) == 1);
        CHECK( fabs( R-1) < 1e-9);
        CHECK( fabs( Z-2) < 1e-14);
    }
}
