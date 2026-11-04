#include <iostream>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "dg/algorithm.h"
#include "magnetic_field.h"
#include "ds.h"
#include "toroidal.h"
#include "catch2/catch_all.hpp"

const double R_0 = 10;
const double a  = 1; //small radius

TEST_CASE( "Fieldaligned boundary in z")
{
    // In a toroidal field the plus and minus transformations just shift planes
    // Independent of the limiter z is not periodic in the communicator
#ifdef WITH_MPI
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {1,1,0}, {false,false,false});
#endif
    const dg::x::CylindricalGrid3d g3d( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        3, 8, 8, 4, dg::NEU, dg::NEU, dg::DIR
#ifdef WITH_MPI
    , comm
#endif
    );
    const double hz = g3d.hz();
    auto func = [=]( double R, double Z, double phi){
        return exp(R-R_0)*(1+Z*Z)*(1+sin(phi));
    };
    const dg::geo::TokamakMagneticField mag = dg::geo::createToroidalField( R_0);
    auto bhat = dg::geo::createBHat(mag);
    const dg::x::DVec fun = dg::evaluate( func, g3d);
    const dg::x::DVec vol3d = dg::create::volume( g3d);
    dg::x::DVec result( fun);
    auto error = [&]( const dg::x::DVec& expected){
        dg::blas1::axpby( 1., expected, -1., result);
        return sqrt( dg::blas2::dot( result, vol3d, result));
    };
    const double eps = 1e-12;
    SECTION( "No limiter wraps planes periodically")
    {
        dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>
            fa( bhat, g3d, dg::NEU, dg::NEU, dg::geo::NoLimiter(), 1e-8, 1, 1, -1, "dg", false);
        fa( dg::geo::einsPlus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            return func( R, Z, phi+hz);}, g3d)) < eps);
        fa( dg::geo::einsMinus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            return func( R, Z, phi-hz);}, g3d)) < eps);
    }
    SECTION( "Full limiter sets ghost values")
    {
        const double left = 1., right = 2.;
        const auto bc = GENERATE( dg::DIR, dg::NEU);
        INFO( "bc "<<dg::bc2str(bc));
        dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>
            fa( bhat, g3d, dg::NEU, dg::NEU, dg::geo::FullLimiter(), 1e-8, 1, 1, -1, "dg", false);
        fa.set_boundaries( bc, left, right);
        fa( dg::geo::einsPlus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi + hz < 2.*M_PI)
                return func( R, Z, phi+hz);
            return bc == dg::DIR ? 2.*right - func( R,Z,phi)
                                 : func( R,Z,phi) + hz*right;
            }, g3d)) < eps);
        fa( dg::geo::einsMinus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi - hz > 0)
                return func( R, Z, phi-hz);
            return bc == dg::DIR ? 2.*left - func( R,Z,phi)
                                 : func( R,Z,phi) - hz*left;
            }, g3d)) < eps);
    }
}

// Provides what the assign_bc functions need from Fieldaligned for a single point
struct SinglePoint
{
    SinglePoint( double delta, double hm, double hp, double bm, double bo, double bp) :
        m_delta(delta), m_hbm(1,hm), m_hbp(1,hp), m_bbm(1,bm), m_bbo(1,bo), m_bbp(1,bp){}
    double deltaPhi() const{ return m_delta;}
    const dg::HVec& hbm() const{ return m_hbm;}
    const dg::HVec& hbp() const{ return m_hbp;}
    const dg::HVec& bbm() const{ return m_bbm;}
    const dg::HVec& bbo() const{ return m_bbo;}
    const dg::HVec& bbp() const{ return m_bbp;}
    private:
    double m_delta;
    dg::HVec m_hbm, m_hbp, m_bbm, m_bbo, m_bbp;
};

TEST_CASE( "Boundary conditions along fieldline")
{
    // s=0 is the center, the boundaries lie at s=-hbm and s=hbp and the minus
    // and plus points at s=-delta and s=+delta (they are outside if the mask says so)
    const double delta = 0.3, hbm = 0.2, hbp = 0.1, eps = 1e-14;
    const double c0 = 1.3, c1 = 0.7, c2 = -1.1;
    auto poly = [=]( double s){ return c0 + c1*s + c2*s*s;};
    auto dpoly = [=]( double s){ return c1 + 2.*c2*s;};
    auto line = [=]( double s){ return c0 + c1*s;};
    enum Mask{ minus, both, plus};
    const Mask mask = GENERATE( minus, both, plus);
    INFO( "Mask "<<(mask == minus ? "minus" : mask == both ? "both" : "plus"));
    const SinglePoint fa( delta, hbm, hbp, mask == minus, mask == both, mask == plus);
    dg::HVec fmg(1), fpg(1);

    SECTION( "2nd order reproduces a parabola")
    {
        const dg::HVec fm(1,poly(-delta)), fo(1,poly(0)), fp(1,poly(delta));
        dg::geo::assign_bc_along_field_2nd( fa, fm, fo, fp, fmg, fpg, dg::DIR,
                {poly(-hbm), poly(hbp)});
        CHECK_THAT( fmg[0], Catch::Matchers::WithinAbs( poly(-delta), eps));
        CHECK_THAT( fpg[0], Catch::Matchers::WithinAbs( poly(+delta), eps));
        dg::geo::assign_bc_along_field_2nd( fa, fm, fo, fp, fmg, fpg, dg::NEU,
                {dpoly(-hbm), dpoly(hbp)});
        CHECK_THAT( fmg[0], Catch::Matchers::WithinAbs( poly(-delta), eps));
        CHECK_THAT( fpg[0], Catch::Matchers::WithinAbs( poly(+delta), eps));
    }
    SECTION( "1st order reproduces a line")
    {
        const dg::HVec fm(1,line(-delta)), fp(1,line(delta));
        dg::geo::assign_bc_along_field_1st( fa, fm, fp, fmg, fpg, dg::DIR,
                {line(-hbm), line(hbp)});
        CHECK_THAT( fmg[0], Catch::Matchers::WithinAbs( line(-delta), eps));
        CHECK_THAT( fpg[0], Catch::Matchers::WithinAbs( line(+delta), eps));
        dg::geo::assign_bc_along_field_1st( fa, fm, fp, fmg, fpg, dg::NEU, {c1, c1});
        CHECK_THAT( fmg[0], Catch::Matchers::WithinAbs( line(-delta), eps));
        CHECK_THAT( fpg[0], Catch::Matchers::WithinAbs( line(+delta), eps));
    }
    SECTION( "Only points outside are replaced")
    {
        // These input values are inconsistent with the boundary values
        const dg::HVec fm(1,5.), fp(1,7.);
        dg::geo::assign_bc_along_field_1st( fa, fm, fp, fmg, fpg, dg::DIR, {1., 2.});
        CHECK( (fmg[0] == 5.) == (mask == plus));
        CHECK( (fpg[0] == 7.) == (mask == minus));
        dg::geo::assign_bc_along_field_1st( fa, fm, fp, fmg, fpg, dg::NEU, {1., 2.});
        CHECK( (fmg[0] == 5.) == (mask != minus));
        CHECK( (fpg[0] == 7.) == (mask != plus));
    }
    SECTION( "Swap perpendicular boundary condition")
    {
        const dg::HVec fm(1,5.), fp(1,7.);
        dg::geo::swap_bc_perp( fa, fm, fp, fmg, fpg);
        CHECK( fmg[0] == (mask == plus ? +5. : -5.));
        CHECK( fpg[0] == (mask == minus ? +7. : -7.));
    }
}
