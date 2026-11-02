#include <iostream>
#include <iomanip>

#include "dg/algorithm.h"

#include "solovev.h"
#include "taylor.h"
#include "magnetic_field.h"
#include "catch2/catch_all.hpp"

struct JPhi
{
    JPhi( dg::geo::solovev::Parameters gp): R_0(gp.R_0), A(gp.A){}
    double operator()(double R, double, double)const
    {
        return ((A-1.)*R - A*R_0*R_0/R)/R_0/R_0/R_0;
    }
    private:
    double R_0, A;
};
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

const std::string geometry_params_taylor = R"asdf({
    "c" :[-0.730004529864617913443309878090712685057289011596696572774289301646,
    -0.013752279867251353918502943858642437963550235876487492144645446786828705,
    -1.043934755061476893780514356717827737689350282272451493765507415754555292,
    0.0652044051606185896009860257256201123783699931775930983650947832207138720,
    -0.075588684108624986707612696927980524120446412329125959454845935222288227,
    -0.166244745975333508158120968358969752878367824897146548390838214529077855,
    0.2339447574961647090247023551091674910845416152319523071184183514964211299,
    -0.106340609874360718259386870009150862837873566733690651435539105985033601,
    -0.195348950338876734547335892506041028973952380901869144068532010742858264,
    0.0185335490893965223839308241407808534672715395435233413969584688870300260,
    3.4373031197519723201231633783929522278128352814841152895441705480656431803,
    4.9016827404530536864875893866884491045861716739824591650253786527604617409],
    "inverseaspectratio" : 0.41071428571428575,
    "elongation"         : 1.75,
    "triangularity"      : 0.25,
    "equilibrium"  : "taylor",
    "description" : "standardX"
})asdf";

TEST_CASE("Solovev description")
{
    // All parameters but "description"
    const std::string params = R"asdf({
        "A" : 0.0, "PP": 1, "PI": 1,
        "c" :[ 1,2,3,4,5,6,7,8,9,10,11,12],
        "R_0" : 10,
        "inverseaspectratio" : 0.25,
        "elongation" : 1,
        "triangularity" : 0
    })asdf";
    const dg::file::error mode = GENERATE( dg::file::error::is_silent,
            dg::file::error::is_warning, dg::file::error::is_throw);
    dg::file::WrappedJsonValue js( dg::file::string2Json( params), mode);
    SECTION( "Non-zero up-down asymmetric coefficients")
    {
        CHECK( dg::geo::solovev::Parameters(js).description == "standardX");
    }
    SECTION( "Up-down symmetric")
    {
        for( unsigned i=7; i<12; i++)
            js.asJson()["c"][i] = 0;
        CHECK( dg::geo::solovev::Parameters(js).description == "standardO");
    }
    SECTION( "Toroidal")
    {
        js.asJson()["PP"] = 0;
        CHECK( dg::geo::solovev::Parameters(js).description == "none");
    }
    SECTION( "Given description is taken")
    {
        js.asJson()["PP"] = 0;
        js.asJson()["description"] = "doubleX";
        CHECK( dg::geo::solovev::Parameters(js).description == "doubleX");
    }
}

TEST_CASE("Taylor")
{
    dg::file::WrappedJsonValue js( dg::file::error::is_silent);
    js.asJson() = dg::file::string2Json( geometry_params_taylor);
    const double R_0 = GENERATE( 1., 10., 547.891714877869);
    js.asJson()["R_0"] = R_0;
    dg::geo::solovev::Parameters gp(js);
    dg::geo::TokamakMagneticField mag = dg::geo::createTaylorField(gp);
    const double R = 1.15*R_0, Z = -0.1*R_0, h = 1e-3*R_0;
    const double lambda = gp.c[11]/R_0;
    INFO( "R_0 "<<R_0);

    // 4th order central differences
    auto dR = [h]( const dg::geo::CylindricalFunctor& f, double r, double z){
        return (-f(r+2*h,z)+8*f(r+h,z)-8*f(r-h,z)+f(r-2*h,z))/(12*h);
    };
    auto dZ = [h]( const dg::geo::CylindricalFunctor& f, double r, double z){
        return (-f(r,z+2*h)+8*f(r,z+h)-8*f(r,z-h)+f(r,z-2*h))/(12*h);
    };
    const double eps = 1e-6;
    SECTION( "Derivatives")
    {
        CHECK_THAT( mag.psipR()(R,Z),  Catch::Matchers::WithinRel( dR( mag.psip(),  R, Z), eps));
        CHECK_THAT( mag.psipZ()(R,Z),  Catch::Matchers::WithinRel( dZ( mag.psip(),  R, Z), eps));
        CHECK_THAT( mag.psipRR()(R,Z), Catch::Matchers::WithinRel( dR( mag.psipR(), R, Z), eps));
        CHECK_THAT( mag.psipRZ()(R,Z), Catch::Matchers::WithinRel( dZ( mag.psipR(), R, Z), eps));
        CHECK_THAT( mag.psipZZ()(R,Z), Catch::Matchers::WithinRel( dZ( mag.psipZ(), R, Z), eps));
        CHECK_THAT( mag.ipolR()(R,Z),  Catch::Matchers::WithinRel( dR( mag.ipol(),  R, Z), eps));
        CHECK_THAT( mag.ipolZ()(R,Z),  Catch::Matchers::WithinRel( dZ( mag.ipol(),  R, Z), eps));
    }
    SECTION( "Grad-Shafranov")
    {
        double gs = mag.psipRR()(R,Z) - mag.psipR()(R,Z)/R + mag.psipZZ()(R,Z);
        CHECK_THAT( gs, Catch::Matchers::WithinRel( -lambda*lambda*mag.psip()(R,Z), eps));
    }
    SECTION( "Force free")
    {
        // Curl in (R,Z,phi) with covariant B_R = B^R, B_Z = B^Z, B_phi = R^2 B^phi
        dg::geo::CylindricalFunctor BR = dg::geo::BFieldR(mag);
        dg::geo::CylindricalFunctor BZ = dg::geo::BFieldZ(mag);
        dg::geo::BFieldP BP(mag);
        dg::geo::CylindricalFunctor Bphi = [&BP]( double r, double z){
            return r*r*BP(r,z);
        };
        CHECK_THAT( dZ( Bphi, R, Z)/R, Catch::Matchers::WithinRel( lambda*BR(R,Z), eps));
        CHECK_THAT( -dR( Bphi, R, Z)/R, Catch::Matchers::WithinRel( lambda*BZ(R,Z), eps));
        CHECK_THAT( (dR( BZ, R, Z) - dZ( BR, R, Z))/R,
                Catch::Matchers::WithinRel( lambda*BP(R,Z), eps));
    }
}

TEST_CASE("Solovev")
{
    dg::file::WrappedJsonValue js( dg::file::error::is_warning);
    js.asJson() = dg::file::string2Json( geometry_params_Xpoint);

    dg::geo::solovev::Parameters gp(js);
    dg::geo::TokamakMagneticField mag = dg::geo::createSolovevField(gp);

    SECTION( "X-point")
    {
        double R_X = gp.R_0-1.1*gp.triangularity*gp.a;
        double Z_X = -1.1*gp.elongation*gp.a;
        REQUIRE( gp.hasXpoint());

        dg::geo::findXpoint( mag.get_psip(), R_X, Z_X);
        INFO(  "X-point found at "<<R_X << " "<<Z_X<<" with Psip "<<mag.psip()(R_X, Z_X));
        INFO("     R - Factor "<<(gp.R_0-R_X)/gp.triangularity/gp.a
              << " Z - factor "<<-(Z_X/gp.elongation/gp.a));
        CHECK( fabs(R_X - (gp.R_0-1.1*gp.triangularity*gp.a) ) < 1e-10);
        CHECK( fabs(Z_X + 1.1*gp.elongation*gp.a ) < 1e-10);
    }

    SECTION( "Accuracy of psi")
    {
        const double R_H = gp.R_0-gp.triangularity*gp.a;
        const double Z_H = gp.elongation*gp.a;
        const double alpha_ = asin(gp.triangularity);
        const double N1 = -(1.+alpha_)/(gp.a*gp.elongation*gp.elongation)*(1.+alpha_);
        const double N2 =  (1.-alpha_)/(gp.a*gp.elongation*gp.elongation)*(1.-alpha_);
        const double N3 = -gp.elongation/(gp.a*cos(alpha_)*cos(alpha_));
        const double R_X = gp.R_0-1.1*gp.triangularity*gp.a;
        const double Z_X = -1.1*gp.elongation*gp.a;

        double result;
        result = mag.psip()(gp.R_0 + gp.a, 0.);
        INFO( "psip( 1+e,0)           "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psip()(gp.R_0 - gp.a, 0.);
        INFO( "psip( 1-e,0)           "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psip()(R_H, Z_H);
        INFO( "psip( 1-de,ke)         "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psip()(R_X, Z_X);
        INFO( "psip( 1-1.1de,-1.1ke)  "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psipZ()(gp.R_0 + gp.a, 0.);
        INFO( "psipZ( 1+e,0)          "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psipZ()(gp.R_0 - gp.a, 0.);
        INFO( "psipZ( 1-e,0)          "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psipR()(R_H,Z_H);
        INFO( "psipR( 1-de,ke)        "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psipR()(R_X,Z_X);
        INFO( "psipR( 1-1.1de,-1.1ke) "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psipZ()(R_X,Z_X);
        INFO( "psipZ( 1-1.1de,-1.1ke) "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psipZZ()(gp.R_0+gp.a,0.)+N1*mag.psipR()(gp.R_0+gp.a,0);
        INFO( "psipZZ( 1+e,0)         "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psipZZ()(gp.R_0-gp.a,0.)+N2*mag.psipR()(gp.R_0-gp.a,0);
        INFO( "psipZZ( 1-e,0)         "<<result);
        CHECK( fabs( result) < 1e-13);
        result = mag.psipRR()(R_H,Z_H)+N3*mag.psipZ()(R_H,Z_H);
        INFO( "psipRR( 1-de,ke)       "<<result);
        CHECK( fabs( result) < 1e-13);
    }

    SECTION( "Accuracy of bhat")
    {
        double Rmin=gp.R_0-gp.a;
        double Zmin=-gp.a*gp.elongation;
        double Rmax=gp.R_0+gp.a;
        double Zmax=gp.a*gp.elongation;
        dg::CylindricalGrid3d grid3d(Rmin,Rmax,Zmin,Zmax, 0, 2.*M_PI, 3,100,100,32);
        dg::geo::CylindricalVectorLvl0 bhat_ = dg::geo::createBHat( mag);
        std::array<dg::HVec, 3> bhat;
        dg::pushForward( bhat_.x(), bhat_.y(), bhat_.z(),
                bhat[0], bhat[1], bhat[2], grid3d);
        std::array<dg::HVec, 3> bhat_covariant(bhat);
        dg::tensor::inv_multiply3d( grid3d.metric(), bhat[0], bhat[1], bhat[2],
                bhat_covariant[0], bhat_covariant[1], bhat_covariant[2]);
        dg::HVec normb( bhat[0]), one3d = dg::evaluate( dg::one, grid3d);
        dg::blas1::pointwiseDot( 1., bhat[0], bhat_covariant[0],
                                 1., bhat[1], bhat_covariant[1],
                                 0., normb);
        dg::blas1::pointwiseDot( 1., bhat[2], bhat_covariant[2],
                                 1., normb);
        dg::blas1::axpby( 1., one3d, -1, normb);
        double error = sqrt(dg::blas1::dot( normb, normb));
        INFO( "Error in norm b == 1 :  "<<error);
        CHECK( error < 1e-12);
    }

    SECTION( "Push Forward curvatures")
    {
        double Rmin=gp.R_0-gp.a;
        double Zmin=-gp.a*gp.elongation;
        double Rmax=gp.R_0+gp.a;
        double Zmax=gp.a*gp.elongation;
        dg::CylindricalGrid3d grid3d(Rmin,Rmax,Zmin,Zmax, 0, 2.*M_PI, 3,100,100,32);
        dg::geo::CylindricalVectorLvl0 curvB_ = dg::geo::createTrueCurvatureNablaB( mag);
        dg::geo::CylindricalVectorLvl0 curvK_ = dg::geo::createTrueCurvatureKappa( mag);
        std::array<dg::HVec, 3> curvB, curvK;
        //Test NablaTimes B = B^2( curvK - curvB)
        dg::pushForward( curvB_.x(), curvB_.y(), curvB_.z(),
                curvB[0], curvB[1], curvB[2], grid3d);
        dg::pushForward( curvK_.x(), curvK_.y(), curvK_.z(),
                curvK[0], curvK[1], curvK[2], grid3d);
        //Test NablaTimes B = B^2( curvK - curvB)
        dg::blas1::axpby( 1., curvK, -1., curvB);
        dg::HVec Bmodule = dg::pullback( dg::geo::Bmodule(mag), grid3d);
        dg::blas1::pointwiseDot( Bmodule, Bmodule, Bmodule);
        for( int i=0; i<3; i++)
            dg::blas1::pointwiseDot( Bmodule, curvB[i], curvB[i]);
        dg::HVec R = dg::pullback( dg::cooX3d, grid3d);
        dg::HVec IR =  dg::pullback( mag.ipolR(), grid3d);
        dg::blas1::pointwiseDivide( gp.R_0, IR, R, 0., IR);
        dg::HVec IZ =  dg::pullback( mag.ipolZ(), grid3d);
        dg::blas1::pointwiseDivide( gp.R_0, IZ, R, 0., IZ);
        dg::HVec IP =  dg::pullback( JPhi( gp), grid3d);
        dg::blas1::pointwiseDivide( gp.R_0, IP, R, 0., IP);
        dg::blas1::axpby( 1., IZ, -1., curvB[0]);
        dg::blas1::axpby( 1., IR, +1., curvB[1]);
        dg::blas1::axpby( 1., IP, -1., curvB[2]);
        for( int i=0; i<3; i++)
        {
            double error = sqrt(dg::blas1::dot( curvB[i], curvB[i] ) );
            INFO( "Error in curv "<<i<<" :   "<<error);
            CHECK( error < 1e-13);
        }
    }
}

