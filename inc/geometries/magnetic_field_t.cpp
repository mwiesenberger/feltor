
#include <iostream>
#include <iomanip>

#include "dg/algorithm.h"

#include "solovev.h"
#include "guenter.h"
#include "polynomial.h"
#include "taylor.h"
#include "toroidal.h"
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

// 4th order central differences
static double diffR( const dg::geo::CylindricalFunctor& f, double R, double Z, double h)
{
    return (-f(R+2*h,Z)+8*f(R+h,Z)-8*f(R-h,Z)+f(R-2*h,Z))/(12*h);
}
static double diffZ( const dg::geo::CylindricalFunctor& f, double R, double Z, double h)
{
    return (-f(R,Z+2*h)+8*f(R,Z+h)-8*f(R,Z-h)+f(R,Z-2*h))/(12*h);
}
static auto within( double value)
{
    return Catch::Matchers::WithinAbs( value, 1e-6*(fabs(value) + 1e-5));
}

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
        const double h = 1e-3;
        dg::geo::CylindricalFunctor bR = dg::geo::ToroidalBHatR(mag);
        dg::geo::CylindricalFunctor bZ = dg::geo::ToroidalBHatZ(mag);
        dg::geo::CylindricalFunctor Rb = [&bR]( double r, double z){ return r*bR(r,z);};
        dg::geo::Btor btor(mag);
        dg::geo::CylindricalFunctor lnB = [&btor]( double r, double z){
            return log( fabs( btor(r,z)));
        };
        const double eps = 1e-6;
        CHECK_THAT( torgradLnB(R,Z), Catch::Matchers::WithinRel(
                bR(R,Z)*diffR( lnB, R, Z, h) + bZ(R,Z)*diffZ( lnB, R, Z, h), eps));
        CHECK_THAT( tordiv(R,Z), Catch::Matchers::WithinRel(
                diffR( Rb, R, Z, h)/R + diffZ( bZ, R, Z, h), eps));
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

// Compare all functors of a magnetic field with finite differences of the
// functor they are derived from
static void check_field( const dg::geo::TokamakMagneticField& mag, double R, double Z, double h)
{
    using namespace dg::geo;
    auto dR = [&]( const CylindricalFunctor& f){ return diffR( f, R, Z, h);};
    auto dZ = [&]( const CylindricalFunctor& f){ return diffZ( f, R, Z, h);};
    // divergence of vector v^R, v^Z (v^phi does not contribute)
    auto div = [&]( const CylindricalFunctor& vR, const CylindricalFunctor& vZ){
        CylindricalFunctor Rv = [&vR]( double r, double z){ return r*vR(r,z);};
        return dR( Rv)/R + dZ( vZ);
    };
    SECTION( "Flux function and current")
    {
        CHECK_THAT( mag.psipR()(R,Z), within( dR( mag.psip())));
        CHECK_THAT( mag.psipZ()(R,Z), within( dZ( mag.psip())));
        CHECK_THAT( mag.psipRR()(R,Z), within( dR( mag.psipR())));
        CHECK_THAT( mag.psipRZ()(R,Z), within( dZ( mag.psipR())));
        CHECK_THAT( mag.psipZZ()(R,Z), within( dZ( mag.psipZ())));
        CHECK_THAT( mag.ipolR()(R,Z), within( dR( mag.ipol())));
        CHECK_THAT( mag.ipolZ()(R,Z), within( dZ( mag.ipol())));
    }
    CylindricalFunctor bR = BHatR(mag), bZ = BHatZ(mag), bP = BHatP(mag);
    SECTION( "Field strength and unit vector")
    {
        CHECK_THAT( BR(mag)(R,Z), within( dR( Bmodule(mag))));
        CHECK_THAT( BZ(mag)(R,Z), within( dZ( Bmodule(mag))));
        CHECK_THAT( ToroidalBR(mag)(R,Z), within( dR( Btor(mag))));
        CHECK_THAT( ToroidalBZ(mag)(R,Z), within( dZ( Btor(mag))));
        CHECK_THAT( BHatRR(mag)(R,Z), within( dR( bR)));
        CHECK_THAT( BHatRZ(mag)(R,Z), within( dZ( bR)));
        CHECK_THAT( BHatZR(mag)(R,Z), within( dR( bZ)));
        CHECK_THAT( BHatZZ(mag)(R,Z), within( dZ( bZ)));
        CHECK_THAT( BHatPR(mag)(R,Z), within( dR( bP)));
        CHECK_THAT( BHatPZ(mag)(R,Z), within( dZ( bP)));
        CHECK_THAT( Divb(mag)(R,Z), within( div( bR, bZ)));
        CHECK_THAT( GradLnB(mag)(R,Z), within( bR(R,Z)*dR( LnB(mag)) + bZ(R,Z)*dZ( LnB(mag))));
        CylindricalFunctor vR = [&]( double r, double z){ return bR(r,z)/bP(r,z);};
        CylindricalFunctor vZ = [&]( double r, double z){ return bZ(r,z)/bP(r,z);};
        CHECK_THAT( DivVVP(mag)(R,Z), within( div( vR, vZ)));
    }
    SECTION( "True curvature")
    {
        // K_kappa = curl(b) / B with (R,Z,phi) right handed and the covariant
        // components b_R = b^R, b_Z = b^Z and b_phi = R^2 b^phi
        CylindricalFunctor bcov = [&bP]( double r, double z){ return r*r*bP(r,z);};
        const double invB = InvB(mag)(R,Z);
        CHECK_THAT( TrueCurvatureKappaR(mag)(R,Z), within( invB*dZ( bcov)/R));
        CHECK_THAT( TrueCurvatureKappaZ(mag)(R,Z), within( -invB*dR( bcov)/R));
        CHECK_THAT( TrueCurvatureKappaP(mag)(R,Z), within( invB*( dR( bZ) - dZ( bR))/R));
        // K_gradB = b x grad(B) / B^2
        const double gradR = dR( Bmodule(mag)), gradZ = dZ( Bmodule(mag));
        CHECK_THAT( TrueCurvatureNablaBR(mag)(R,Z), within( -R*bP(R,Z)*gradZ*invB*invB));
        CHECK_THAT( TrueCurvatureNablaBZ(mag)(R,Z), within( +R*bP(R,Z)*gradR*invB*invB));
        CHECK_THAT( TrueCurvatureNablaBP(mag)(R,Z), within(
            (bR(R,Z)*gradZ - bZ(R,Z)*gradR)*invB*invB/R));
        CHECK_THAT( TrueDivCurvatureKappa(mag)(R,Z), within(
            div( TrueCurvatureKappaR(mag), TrueCurvatureKappaZ(mag))));
        CHECK_THAT( TrueDivCurvatureNablaB(mag)(R,Z), within(
            div( TrueCurvatureNablaBR(mag), TrueCurvatureNablaBZ(mag))));
    }
    SECTION( "Approximate curvature")
    {
        const int sign = GENERATE( +1, -1);
        CHECK_THAT( DivCurvatureKappa(mag,sign)(R,Z), within(
            div( CurvatureKappaR(mag,sign), CurvatureKappaZ(mag,sign))));
        CHECK_THAT( DivCurvatureNablaB(mag,sign)(R,Z), within(
            div( CurvatureNablaBR(mag,sign), CurvatureNablaBZ(mag,sign))));
        CHECK_THAT( ToroidalDivCurvatureKappa(mag)(R,Z), within(
            div( ToroidalCurvatureKappaR(mag), ToroidalCurvatureKappaZ(mag))));
    }
}

TEST_CASE("Magnetic field functors against finite differences")
{
    SECTION( "Solovev")
    {
        dg::file::WrappedJsonValue js( dg::file::string2Json( geometry_params_Xpoint),
                dg::file::error::is_silent);
        js.asJson()["A"] = 1.;
        dg::geo::solovev::Parameters gp(js);
        check_field( dg::geo::createSolovevField(gp), 0.615, -0.0525, 1e-3);
    }
    SECTION( "Polynomial")
    {
        dg::file::WrappedJsonValue js( dg::file::string2Json( R"asdf({
            "PP": -1.5, "PI": 2, "R_0": 4, "M": 4, "N": 4,
            "c": [ 1.0, -0.5,  0.3,  0.2,
                  -0.7,  0.4,  0.6, -0.1,
                   0.8,  0.9, -0.3,  0.5,
                  -0.2,  0.7,  0.1, -0.4]
        })asdf"), dg::file::error::is_silent);
        dg::geo::polynomial::Parameters gp(js);
        check_field( dg::geo::createPolynomialField(gp), 4.3, 0.4, 1e-3);
    }
    SECTION( "Taylor")
    {
        dg::file::WrappedJsonValue js( dg::file::string2Json( R"asdf({
            "R_0": 10,
            "c": [ -0.7300045298646179, -0.01375227986725135, -1.0439347550614769,
                0.06520440516061859, -0.07558868410862499, -0.1662447459753335,
                0.2339447574961647, -0.10634060987436072, -0.19534895033887673,
                0.01853354908939652, 3.437303119751972, 4.901682740453054]
        })asdf"), dg::file::error::is_silent);
        dg::geo::solovev::Parameters gp(js);
        check_field( dg::geo::createTaylorField(gp), 11.5, -1.0, 1e-3);
    }
    SECTION( "Guenter")
    {
        check_field( dg::geo::createGuenterField( 10, 20), 10.3, 0.4, 1e-3);
    }
    SECTION( "Circular")
    {
        check_field( dg::geo::createCircularField( 10, 20, 1.2, 0.8), 10.3, 0.4, 1e-3);
    }
}

TEST_CASE("Periodify")
{
    // Mirror (3.5,0.3) at R=4 or shift (2.5,0.3) by the box length
    auto mag = dg::geo::createCircularField( 3, 10, 1, 1);
    const double R = 4.5, Z = 0.3;
    const double psi = 1 - 0.25 - 0.09;
    SECTION( "NEU")
    {
        auto per = dg::geo::periodify( mag, 2, 4, -1, 1, dg::NEU, dg::NEU);
        CHECK( fabs( per.psip()(R,Z) - psi) < 1e-14);
        CHECK( fabs( per.ipol()(R,Z) - 10) < 1e-14);
        CHECK_THAT( per.psipR()(R,Z), within( diffR( per.psip(), R, Z, 1e-3)));
        CHECK_THAT( per.psipZ()(R,Z), within( diffZ( per.psip(), R, Z, 1e-3)));
        CHECK_THAT( per.psipRR()(R,Z), within( diffR( per.psipR(), R, Z, 1e-3)));
    }
    SECTION( "DIR")
    {
        auto per = dg::geo::periodify( mag, 2, 4, -1, 1, dg::DIR, dg::NEU);
        CHECK( fabs( per.psip()(R,Z) + psi) < 1e-14);
        CHECK_THAT( per.psipR()(R,Z), within( diffR( per.psip(), R, Z, 1e-3)));
        CHECK_THAT( per.psipRR()(R,Z), within( diffR( per.psipR(), R, Z, 1e-3)));
    }
    SECTION( "PER")
    {
        auto per = dg::geo::periodify( mag, 2, 4, -1, 1, dg::PER, dg::PER);
        CHECK( fabs( per.psip()(R,Z) - psi) < 1e-14);
        CHECK_THAT( per.psipR()(R,Z), within( diffR( per.psip(), R, Z, 1e-3)));
    }
}

TEST_CASE("Alignment and projection tensors")
{
    dg::file::WrappedJsonValue js( dg::file::string2Json( geometry_params_Xpoint),
            dg::file::error::is_silent);
    js.asJson()["A"] = 1.;
    dg::geo::solovev::Parameters gp(js);
    dg::geo::TokamakMagneticField mag = dg::geo::createSolovevField(gp);
    dg::CylindricalGrid3d grid( gp.R_0-gp.a, gp.R_0+gp.a, -gp.a, gp.a, 0, 2.*M_PI,
            3, 10, 10, 4);
    auto bhat = dg::geo::createBHat( mag);
    std::array<dg::HVec,3> b, bcov, result;
    dg::pushForward( bhat.x(), bhat.y(), bhat.z(), b[0], b[1], b[2], grid);
    bcov = result = b;
    dg::tensor::inv_multiply3d( grid.metric(), b[0], b[1], b[2], bcov[0], bcov[1], bcov[2]);
    SECTION( "Alignment tensor leaves b unchanged")
    {
        auto chi = dg::geo::createAlignmentTensor( bhat, grid);
        dg::tensor::multiply3d( chi, bcov[0], bcov[1], bcov[2], result[0], result[1], result[2]);
        for( unsigned u=0; u<3; u++)
        {
            dg::blas1::axpby( 1., b[u], -1., result[u]);
            INFO( "Component "<<u);
            CHECK( sqrt( dg::blas1::dot( result[u], result[u])) < 1e-12);
        }
    }
    SECTION( "Projection tensor removes b")
    {
        auto h = dg::geo::createProjectionTensor( bhat, grid);
        dg::tensor::multiply3d( h, bcov[0], bcov[1], bcov[2], result[0], result[1], result[2]);
        for( unsigned u=0; u<3; u++)
        {
            INFO( "Component "<<u);
            CHECK( sqrt( dg::blas1::dot( result[u], result[u])) < 1e-12);
        }
    }
}
