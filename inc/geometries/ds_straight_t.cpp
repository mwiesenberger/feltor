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
