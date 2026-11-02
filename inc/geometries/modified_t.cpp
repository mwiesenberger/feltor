#include <iostream>


#include "dg/file/json_utilities.h"
#include "make_field.h"
#include "modified.h"
#include "catch2/catch_all.hpp"

/* Input-File for COMPASS axisymmetric solovev equilibrium *
  I const with X-point, Te_50eV,B0_1T,deuterium, Cref~2e-5
  ----------------------------------------------------------- */

const std::string geometry_params_toroidal = R"asdf({
    "R_0": 100,
    "equilibrium" : "toroidal"
})asdf";

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
const std::string geometry_params_2Xpoint = R"asdf({
    "A": 0,
    "PP": 1,
    "PI": 1,
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
    "R_0" : 545,
    "inverseaspectratio": 0.3211009174311926,
    "elongation": 1.44,
    "triangularity": 0.3,
    "equilibrium": "solovev",
    "description" : "doubleX"

})asdf";

const std::string geometry_params_circular = R"asdf({
    "equilibrium" : "circular",
    "I_0" : 10,
    "R_0" : 3,
    "a" : 1,
    "b" : 1
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
    "R_0"                : 547.891714877869,
    "inverseaspectratio" : 0.41071428571428575,
    "elongation"         : 1.75,
    "triangularity"      : 0.25,
    "equilibrium"  : "taylor",
    "description" : "standardX"
})asdf";
const std::string geometry_params_polynomial = R"asdf({
    "PP": 1,
    "PI": 1,
    "c" :[  1,2,3,4],
    "R_0"                : 4,
    "M" : 2,
    "N" : 2,
    "equilibrium"  : "polynomial",
    "description" : "standardX"
})asdf";
const std::string geometry_params_guenter = R"asdf({
    "equilibrium" : "guenter",
    "I_0" : 10,
    "R_0" : 3
})asdf";

// const std::array<double,2> xpoint = {431.552779491927, -433.176887075315};

//@ ----------------------------------------------------------

static dg::file::WrappedJsonValue solovev_params( const std::string& params, double PP)
{
    dg::file::WrappedJsonValue js( dg::file::string2Json( params),
            dg::file::error::is_silent);
    js.asJson()["PP"] = PP;
    return js;
}

// The wall must be 1 where the modified flux function and the current are
// constant and 0 where they are not modified
static void check_wall( const dg::file::WrappedJsonValue& js,
    const std::string& modifier, double R0, double R1, double Z0, double Z1)
{
    dg::file::WrappedJsonValue jsmod( dg::file::string2Json( modifier),
            dg::file::error::is_silent);
    dg::geo::CylindricalFunctor wall, transition;
    dg::geo::TokamakMagneticField mod = dg::geo::createModifiedField( js, jsmod,
            wall, transition);
    dg::geo::TokamakMagneticField mag = dg::geo::createMagneticField( js);
    dg::Grid2d grid( R0, R1, Z0, Z1, 1, 100, 100);
    const dg::HVec chi = dg::evaluate( wall, grid);
    const dg::HVec modR = dg::evaluate( mod.psipR(), grid);
    const dg::HVec modZ = dg::evaluate( mod.psipZ(), grid);
    const dg::HVec magR = dg::evaluate( mag.psipR(), grid);
    const dg::HVec magZ = dg::evaluate( mag.psipZ(), grid);
    const dg::HVec modIR = dg::evaluate( mod.ipolR(), grid);
    const dg::HVec modIZ = dg::evaluate( mod.ipolZ(), grid);
    const dg::HVec magIR = dg::evaluate( mag.ipolR(), grid);
    const dg::HVec magIZ = dg::evaluate( mag.ipolZ(), grid);
    // At the edge of the transition the wall rounds to exactly 1 or 0 before the
    // derivatives do
    const double eps = 1e-14;
    unsigned in_wall = 0, outside = 0, in_wall_wrong = 0, outside_wrong = 0;
    for( unsigned i=0; i<chi.size(); i++)
    {
        if( chi[i] == 1)
        {
            in_wall++;
            if( fabs( modR[i]) > eps || fabs( modZ[i]) > eps
                    || fabs( modIR[i]) > eps || fabs( modIZ[i]) > eps)
                in_wall_wrong++;
        }
        else if( chi[i] == 0)
        {
            outside++;
            if( fabs( modR[i] - magR[i]) > eps || fabs( modZ[i] - magZ[i]) > eps
                    || fabs( modIR[i] - magIR[i]) > eps
                    || fabs( modIZ[i] - magIZ[i]) > eps)
                outside_wrong++;
        }
    }
    INFO( modifier);
    CHECK( in_wall > 0);
    CHECK( outside > 0);
    CHECK( in_wall_wrong == 0);
    CHECK( outside_wrong == 0);
}

TEST_CASE( "Equilibrium of magnetic field")
{
    auto get_equilibrium = []( const std::string& params){
        dg::file::WrappedJsonValue js( dg::file::string2Json( params),
                dg::file::error::is_silent);
        return dg::geo::createMagneticField( js).params().getEquilibrium();
    };
    using dg::geo::equilibrium;
    CHECK( get_equilibrium( geometry_params_Xpoint) == equilibrium::solovev);
    CHECK( get_equilibrium( geometry_params_taylor) == equilibrium::taylor);
    CHECK( get_equilibrium( geometry_params_polynomial) == equilibrium::polynomial);
    CHECK( get_equilibrium( geometry_params_guenter) == equilibrium::guenter);
    CHECK( get_equilibrium( geometry_params_toroidal) == equilibrium::toroidal);
    CHECK( get_equilibrium( geometry_params_circular) == equilibrium::circular);
}

// Without modification the magnetic field must stay the same
TEST_CASE( "Modifier none keeps the magnetic field")
{
    const std::string params = GENERATE( geometry_params_Xpoint,
            geometry_params_taylor, geometry_params_circular);
    dg::file::WrappedJsonValue js( dg::file::string2Json( params),
            dg::file::error::is_silent);
    dg::file::WrappedJsonValue jsmod( dg::file::string2Json(
            R"({"type":"none"})"), dg::file::error::is_silent);
    dg::geo::CylindricalFunctor wall, transition;
    dg::geo::TokamakMagneticField mod = dg::geo::createModifiedField( js, jsmod,
            wall, transition);
    dg::geo::TokamakMagneticField mag = dg::geo::createMagneticField( js);
    const double R = 1.15*mag.R0(), Z = 0.1*mag.R0();
    INFO( params);
    CHECK( mod.params().getEquilibrium() == mag.params().getEquilibrium());
    CHECK( mod.psip()(R,Z) == mag.psip()(R,Z));
    CHECK( mod.psipR()(R,Z) == mag.psipR()(R,Z));
    CHECK( mod.psipZ()(R,Z) == mag.psipZ()(R,Z));
    CHECK( mod.psipRR()(R,Z) == mag.psipRR()(R,Z));
    CHECK( mod.psipRZ()(R,Z) == mag.psipRZ()(R,Z));
    CHECK( mod.psipZZ()(R,Z) == mag.psipZZ()(R,Z));
    CHECK( mod.ipol()(R,Z) == mag.ipol()(R,Z));
    CHECK( mod.ipolR()(R,Z) == mag.ipolR()(R,Z));
    CHECK( mod.ipolZ()(R,Z) == mag.ipolZ()(R,Z));
}

TEST_CASE( "Wall region of modified field")
{
    const dg::file::WrappedJsonValue circular( dg::file::string2Json(
            geometry_params_circular), dg::file::error::is_silent);
    const double PP = GENERATE( +1., -1.);
    INFO( "PP "<<PP);
    const auto Xpoint  = solovev_params( geometry_params_Xpoint, PP);
    const auto Xpoint2 = solovev_params( geometry_params_2Xpoint, PP);
    const dg::file::WrappedJsonValue taylor( dg::file::string2Json(
            geometry_params_taylor), dg::file::error::is_silent);
    SECTION( "heaviside")
    {
        check_wall( Xpoint, R"({"type":"heaviside","boundary":1.1,"alpha":0.2})",
                300, 800, -500, 500);
        check_wall( taylor, R"({"type":"heaviside","boundary":1.1,"alpha":0.2})",
                300, 800, -500, 500);
        check_wall( circular, R"({"type":"heaviside","boundary":0.8,"alpha":0.1})",
                2, 4, -1, 1);
    }
    SECTION( "sol_pfr")
    {
        const std::string mod = R"({"type":"sol_pfr","boundary":[1.1,0.97],"alpha":[0.1,0.1]})";
        check_wall( Xpoint, mod, 300, 800, -500, 500);
        check_wall( Xpoint2, mod, 300, 800, -500, 500);
        check_wall( taylor, mod, 300, 800, -500, 500);
    }
    SECTION( "sol_pfr without X-point")
    {
        // the boundaries are values of psi here
        check_wall( circular,
            R"({"type":"sol_pfr","boundary":[0.8,0.2],"alpha":[0.05,0.05]})",
            2, 4, -1, 1);
    }
    SECTION( "sol_pfr_2X")
    {
        check_wall( Xpoint2, R"({"type":"sol_pfr_2X",
            "boundary":[1.1,0.998,1.02,0.97],"alpha":[0.1,0.1,0.05,0.05]})",
            300, 800, -500, 500);
    }
}

TEST_CASE( "Modified flux function")
{
    dg::file::WrappedJsonValue js( dg::file::error::is_warning);
    js.asJson() = dg::file::string2Json( geometry_params_Xpoint);
    dg::geo::TokamakMagneticField mag = dg::geo::createMagneticField( js);
    // 4th order central differences
    const double h = 0.1;
    auto dR = [h]( const dg::geo::CylindricalFunctor& f, double R, double Z){
        return (-f(R+2*h,Z)+8*f(R+h,Z)-8*f(R-h,Z)+f(R-2*h,Z))/(12*h);
    };
    auto dZ = [h]( const dg::geo::CylindricalFunctor& f, double R, double Z){
        return (-f(R,Z+2*h)+8*f(R,Z+h)-8*f(R,Z-h)+f(R,Z-2*h))/(12*h);
    };
    const double R = 650, Z = -100;
    const double psip = mag.psip()(R,Z);
    const double alpha = 0.2*fabs(psip), psi0 = psip + 0.3*alpha;
    const double sign = GENERATE( +1., -1.);
    INFO( "sign "<<sign);

    SECTION( "Derivatives in transition region")
    {
        auto mod = dg::geo::mod::createPsip( dg::geo::mod::everywhere,
                mag.get_psip(), psi0, alpha, sign);
        CHECK( mod.f()(R,Z) != psip); // is modified
        const double eps = 1e-6;
        CHECK_THAT( mod.dfx()(R,Z),  Catch::Matchers::WithinRel( dR( mod.f(),   R, Z), eps));
        CHECK_THAT( mod.dfy()(R,Z),  Catch::Matchers::WithinRel( dZ( mod.f(),   R, Z), eps));
        CHECK_THAT( mod.dfxx()(R,Z), Catch::Matchers::WithinRel( dR( mod.dfx(), R, Z), eps));
        CHECK_THAT( mod.dfxy()(R,Z), Catch::Matchers::WithinRel( dZ( mod.dfx(), R, Z), eps));
        CHECK_THAT( mod.dfxy()(R,Z), Catch::Matchers::WithinRel( dR( mod.dfy(), R, Z), eps));
        CHECK_THAT( mod.dfyy()(R,Z), Catch::Matchers::WithinRel( dZ( mod.dfy(), R, Z), eps));
    }
    SECTION( "Magnetic transition is non-negative for both signs")
    {
        dg::geo::mod::MagneticTransition transition( dg::geo::mod::everywhere,
                mag.psip(), psi0, alpha, sign);
        CHECK( transition(R,Z) > 0);
    }
    SECTION( "Predicate is false")
    {
        auto mod = dg::geo::mod::createPsip( dg::geo::mod::nowhere,
                mag.get_psip(), psi0, alpha, sign);
        CHECK( mod.f()(R,Z) == psip);
        CHECK( mod.dfx()(R,Z) == mag.psipR()(R,Z));
        CHECK( mod.dfy()(R,Z) == mag.psipZ()(R,Z));
        CHECK( mod.dfxx()(R,Z) == mag.psipRR()(R,Z));
        CHECK( mod.dfxy()(R,Z) == mag.psipRZ()(R,Z));
        CHECK( mod.dfyy()(R,Z) == mag.psipZZ()(R,Z));
    }
}

TEST_CASE( "Modified")
{
    SECTION( "Predicate modifier functions")
    {
        CHECK( dg::geo::mod::everywhere( 1,2) == true);
        CHECK( dg::geo::mod::nowhere( 1,2) == false);
        dg::geo::mod::HeavisideZ heavy( 0.8, +1);
        CHECK( heavy( 1,2) == true);
        CHECK( heavy( 1,0.79) == false);
        dg::geo::mod::HeavisideZ heavyM( 0.8, -1);
        CHECK( heavyM( 1,2) == false);
        CHECK( heavyM( 1,0.79) == true);

        dg::geo::mod::RightSideOf twoP( {-1,-2}, {3,0});
        CHECK( twoP( 0,0) == false);
        CHECK( twoP( 0,3) == false);
        CHECK( twoP( 0,-2) == true);

        dg::geo::mod::RightSideOf threeP( {-1,-2}, {3,0}, {0,3});
        CHECK( threeP( 0,0) == false);
        CHECK( threeP( 0,4) == true);
        CHECK( threeP( 0,-2) == true);

        dg::geo::mod::RightSideOf twoM( {3,0}, {-1,-2});
        CHECK( twoM( 0,0) == true);
        CHECK( twoM( 0,3) == true);
        CHECK( twoM( 0,-2) == false);

        dg::geo::mod::RightSideOf threeM( {0,3}, {3,0}, {-1,-2});
        CHECK( threeM( 0,0) == true);
        CHECK( threeM( 0,4) == false);
        CHECK( threeM( 0,-2) == false);

        dg::geo::mod::Above above( {0,0.8}, {0,1.8}); //  same as HeavisideZ
        CHECK( above( 1,2) == true);
        CHECK( above( 1,0.79) == false);
        dg::geo::mod::Above aboveM( {0,0.8}, {0,1.8}, false); //  same as HeavisideZ
        CHECK( aboveM( 1,2) == false);
        CHECK( aboveM( 1,0.79) == true);
    }

    SECTION( "ClosedFieldline region toroidal")
    {
        dg::file::WrappedJsonValue js( dg::file::error::is_warning);
        js.asJson() = dg::file::string2Json( geometry_params_toroidal);
        dg::geo::TokamakMagneticField mag = dg::geo::createMagneticField(
                js);
        dg::geo::mod::ClosedFieldlineRegion closed(mag);
        CHECK( not closed( 100, 0));
        CHECK( not closed( 100, 10));
        dg::geo::mod::ClosedFieldlineRegion open(mag, false);
        CHECK( open( 100, 0));
        CHECK( open( 100, 10));

    }

    SECTION( "ClosedFieldline region 1 X-point")
    {
        dg::file::WrappedJsonValue js( dg::file::error::is_warning);
        js.asJson() = dg::file::string2Json( geometry_params_Xpoint);
        dg::geo::TokamakMagneticField neg_mag = dg::geo::createMagneticField(
                js);
        CHECK( neg_mag.psip()( 550, 0) < 0);
        js.asJson()["PP"] = -1;
        dg::geo::TokamakMagneticField pos_mag = dg::geo::createMagneticField(
                js);
        CHECK( pos_mag.psip()( 550, 0) > 0);
        for( auto mag : {neg_mag, pos_mag})
        {
            dg::geo::mod::ClosedFieldlineRegion closed(mag);
            CHECK(     closed( 550, 0)); // O-point
            CHECK(     closed( 450, -400)); // above X-point
            CHECK( not closed( 350, -360)); // SOL
            CHECK( not closed( 426, -460)); // PFR
            CHECK( not closed( 600, -360)); // SOL
            CHECK( not closed( 700,  360)); // SOL
            dg::geo::mod::ClosedFieldlineRegion open(mag, false);
            CHECK( not open( 550, 0)); // O-point
            CHECK( not open( 450, -400)); // above X-point
            CHECK(     open( 350, -360)); // SOL
            CHECK(     open( 426, -460)); // PFR
            CHECK(     open( 600, -360)); // SOL
            CHECK(     open( 700,  360)); // SOL
        }
    }
    SECTION( "SOL region")
    {
        dg::file::WrappedJsonValue js( dg::file::error::is_warning);
        js.asJson() = dg::file::string2Json( geometry_params_Xpoint);
        dg::geo::TokamakMagneticField mag = dg::geo::createMagneticField( js);
        const dg::geo::mod::SOLRegion sol( mag, dg::geo::Constant(0));
        CHECK( not sol( 550, 0)); // O-point
        CHECK( not sol( 450, -400)); // above X-point
        CHECK( sol( 350, -360)); // SOL
        CHECK( sol( 426, -460)); // PFR
        const dg::geo::mod::SOLRegion wall( mag, dg::geo::Constant(1));
        CHECK( not wall( 350, -360));
        CHECK( not wall( 550, 0));
    }
    SECTION("ClosedFieldline region 2 X-point")
    {
        dg::file::WrappedJsonValue js( dg::file::error::is_warning);
        js.asJson() = dg::file::string2Json( geometry_params_2Xpoint);
        dg::geo::TokamakMagneticField neg_mag = dg::geo::createMagneticField(
                js);
        CHECK( neg_mag.psip()( 550, 0) < 0);
        js.asJson()["PP"] = -1;
        dg::geo::TokamakMagneticField pos_mag = dg::geo::createMagneticField(
                js);
        CHECK( pos_mag.psip()( 550, 0) > 0);
        for( auto mag : {neg_mag, pos_mag})
        {
            dg::geo::mod::ClosedFieldlineRegion closed(mag);
            CHECK(     closed( 550, 0)); // O-point
            CHECK(     closed( 460, -320)); // above 1st X-point
            CHECK(     closed( 460, +320)); // below 2nd X-point
            CHECK( not closed( 350, -360)); // SOL
            CHECK( not closed( 450, -400)); // PFR 1
            CHECK( not closed( 450, +400)); // PFR 2
            CHECK( not closed( 600, -360)); // SOL
            CHECK( not closed( 700,  360)); // SOL
            dg::geo::mod::ClosedFieldlineRegion open(mag, false);
            CHECK( not open( 550, 0)); // O-point
            CHECK( not open( 460, -320)); // above 1st X-point
            CHECK( not open( 460, +320)); // below 2nd X-point
            CHECK(     open( 350, -360)); // SOL
            CHECK(     open( 450, -400)); // PFR 1
            CHECK(     open( 450, +400)); // PFR 2
            CHECK(     open( 600, -360)); // SOL
            CHECK(     open( 700,  360)); // SOL
        }
    }
}
