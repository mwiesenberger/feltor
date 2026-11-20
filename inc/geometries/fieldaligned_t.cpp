#include <iostream>

#ifdef WITH_MPI
#include <mpi.h>
#include <cstring>
#include <cmath>
#endif

#include "dg/algorithm.h"
#include "magnetic_field.h"
#include "fieldaligned.h"
#ifdef WITH_MPI
#include "mpi_fieldaligned.h"
#endif
#include "toroidal.h"
#include "catch2/catch_all.hpp"

const double R_0 = 10;
const double a  = 1; //small radius

TEST_CASE( "Fieldaligned boundary in z")
{
    // In a toroidal field the plus and minus transformations just shift planes
    // Independent of the limiter z is not periodic in the communicator
    // 8 planes so that every process has one with up to 8 processes in z
#ifdef WITH_MPI
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {1,1,0}, {false,false,false});
#endif
    const dg::x::CylindricalGrid3d g3d( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        3, 8, 8, 8, dg::NEU, dg::NEU, dg::DIR
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
    SECTION( "Boundary values from 2d vectors")
    {
        auto left2d = [=]( double R, double Z){ return 1. + R + Z;};
        auto right2d = [=]( double R, double Z){ return 2. + R*Z;};
        dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>
            fa( bhat, g3d, dg::NEU, dg::NEU, dg::geo::FullLimiter(), 1e-8, 1, 1, -1, "dg", false);
        const dg::ClonePtr<dg::x::aGeometry2d> perp_ptr = g3d.perp_grid();
        const dg::x::DVec left = dg::evaluate( left2d, *perp_ptr);
        const dg::x::DVec right = dg::evaluate( right2d, *perp_ptr);
        fa.set_boundaries( dg::DIR, left, right);
        fa( dg::geo::einsPlus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi + hz < 2.*M_PI)
                return func( R, Z, phi+hz);
            return 2.*right2d(R,Z) - func( R,Z,phi);
            }, g3d)) < eps);
        fa( dg::geo::einsMinus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi - hz > 0)
                return func( R, Z, phi-hz);
            return 2.*left2d(R,Z) - func( R,Z,phi);
            }, g3d)) < eps);
    }
    SECTION( "Boundary values from a 3d vector")
    {
        const double scal_left = 0.5, scal_right = 3.;
        dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>
            fa( bhat, g3d, dg::NEU, dg::NEU, dg::geo::FullLimiter(), 1e-8, 1, 1, -1, "dg", false);
        fa.set_boundaries( dg::DIR, fun, scal_left, scal_right);
        fa( dg::geo::einsPlus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi + hz < 2.*M_PI)
                return func( R, Z, phi+hz);
            return (2.*scal_right - 1.)*func( R,Z,phi);
            }, g3d)) < eps);
        fa( dg::geo::einsMinus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi - hz > 0)
                return func( R, Z, phi-hz);
            return (2.*scal_left - 1.)*func( R,Z,phi);
            }, g3d)) < eps);
    }
    SECTION( "Limiter that is not 1 everywhere")
    {
        // Ghost values only where the limiter is 1, periodic elsewhere
        auto limiter = [=]( double R, double){ return R > R_0 ? 1. : 0.;};
        const double right = 2.;
        dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>
            fa( bhat, g3d, dg::NEU, dg::NEU, limiter, 1e-8, 1, 1, -1, "dg", false);
        fa.set_boundaries( dg::DIR, 0., right);
        fa( dg::geo::einsPlus, fun, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi + hz < 2.*M_PI || limiter(R,Z) == 0)
                return func( R, Z, phi+hz);
            return 2.*right - func( R,Z,phi);
            }, g3d)) < eps);
    }
    SECTION( "Ghost values do not depend on the interpolation method")
    {
        const std::string method = GENERATE( "dg", "linear", "cubic",
            "linear-nearest", "cubic-nearest", "dg-linear", "linear-linear");
        const auto bc = GENERATE( dg::DIR, dg::NEU);
        INFO( "method "<<method<<" bc "<<dg::bc2str(bc));
        const double left = 1., right = 2.;
        dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>
            fa( bhat, g3d, dg::NEU, dg::NEU, dg::geo::FullLimiter(), 1e-8, 3, 3, -1, method, false);
        fa.set_boundaries( bc, left, right);
        fa( dg::geo::einsPlus, fun, result);
        // The interpolation in the other planes is not exact
        const dg::x::DVec last = dg::evaluate( [=]( double, double, double phi){
            return phi + hz < 2.*M_PI ? 0. : 1.;}, g3d);
        dg::blas1::pointwiseDot( last, result, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi + hz < 2.*M_PI)
                return 0.;
            return bc == dg::DIR ? 2.*right - func( R,Z,phi)
                                 : func( R,Z,phi) + hz*right;
            }, g3d)) < eps);
        fa( dg::geo::einsMinus, fun, result);
        const dg::x::DVec first = dg::evaluate( [=]( double, double, double phi){
            return phi - hz > 0 ? 0. : 1.;}, g3d);
        dg::blas1::pointwiseDot( first, result, result);
        CHECK( error( dg::evaluate( [=]( double R, double Z, double phi){
            if( phi - hz > 0)
                return 0.;
            return bc == dg::DIR ? 2.*left - func( R,Z,phi)
                                 : func( R,Z,phi) - hz*left;
            }, g3d)) < eps);
    }
}

#ifdef WITH_MPI
TEST_CASE( "Fieldaligned needs a plane on every process")
{
    int size;
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    if( size == 1) // one process always owns all planes
        return;
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {1,1,0}, {false,false,false});
    // one plane fewer than processes in z
    const dg::x::CylindricalGrid3d g3d( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        3, 8, 8, size-1, dg::NEU, dg::NEU, dg::DIR, comm);
    auto bhat = dg::geo::createBHat( dg::geo::createToroidalField( R_0));
    using FA = dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>;
    CHECK_THROWS_AS( FA( bhat, g3d, dg::NEU, dg::NEU, dg::geo::NoLimiter(),
        1e-8, 1, 1, -1, "dg", false), dg::Error);
}
#endif

TEST_CASE( "Fieldaligned distance to boundary")
{
    // In the field v = (0,c,1) fieldlines are straight with dZ/dphi = c
#ifdef WITH_MPI
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {false,false,true});
#endif
    const dg::x::CylindricalGrid3d g3d( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        3, 8, 8, 4, dg::NEU, dg::NEU, dg::PER
#ifdef WITH_MPI
    , comm
#endif
    );
    const double c = 0.5, delta = g3d.hz();
    dg::geo::CylindricalVectorLvl1 vec( dg::geo::Constant(0), dg::geo::Constant(c),
        dg::geo::Constant(1), dg::geo::Constant(0), dg::geo::Constant(0));
    dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>
        fa( vec, g3d, dg::NEU, dg::NEU, dg::geo::NoLimiter(), 1e-10, 1, 1, -1, "dg", false);
    auto outP = [=]( double Z){ return Z + c*delta > a;};
    auto outM = [=]( double Z){ return Z - c*delta < -a;};
    const dg::x::DVec vol = dg::create::volume( g3d);
    auto error = [&]( const dg::x::DVec& value, auto expected){
        dg::x::DVec diff = dg::evaluate( [=]( double, double Z, double){
            return expected(Z);}, g3d);
        dg::blas1::axpby( 1., value, -1., diff);
        return sqrt( dg::blas2::dot( diff, vol, diff));
    };
    REQUIRE( c*delta < a); // no point leaves in both directions
    CHECK( error( fa.hbp(), [=]( double Z){ return outP(Z) ? (a-Z)/c : delta;}) < 1e-6);
    CHECK( error( fa.hbm(), [=]( double Z){ return outM(Z) ? (Z+a)/c : delta;}) < 1e-6);
    CHECK( error( fa.bbp(), [=]( double Z){ return outP(Z) ? 1. : 0.;}) == 0);
    CHECK( error( fa.bbm(), [=]( double Z){ return outM(Z) ? 1. : 0.;}) == 0);
    CHECK( error( fa.bbo(), [=]( double){ return 0.;}) == 0);
    CHECK( error( fa.bphi(), [=]( double){ return 1.;}) == 0);
    CHECK( fa.deltaPhi() == delta);
}

TEST_CASE( "Fieldaligned evaluate")
{
    // In a toroidal field the plus and minus transformations just shift planes
    // so all that remains is the modulation
#ifdef WITH_MPI
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {1,1,0}, {false,false,true});
#endif
    // With MPI the number of planes is not necessarily divisible by the number of ranks
    // (but every rank needs a plane: at least 8 planes for up to 8 ranks)
    const unsigned Nz = GENERATE( 9u, 11u, 13u);
    INFO( "Nz "<<Nz);
    const dg::x::CylindricalGrid3d g3d( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        3, 8, 8, Nz, dg::NEU, dg::NEU, dg::PER
#ifdef WITH_MPI
    , comm
#endif
    );
    const dg::geo::TokamakMagneticField mag = dg::geo::createToroidalField( R_0);
    auto bhat = dg::geo::createBHat(mag);
    dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>
        fa( bhat, g3d, dg::NEU, dg::NEU, dg::geo::NoLimiter(), 1e-8, 1, 1, -1, "dg", false);
    const double hz = g3d.hz();
    const unsigned p0 = 1;
    auto binary = [=]( double R, double Z){ return exp(R-R_0)*(1+Z*Z);};
    // The value at 0 is not 1
    auto unary = [=]( double phi){ return 2.*exp( -phi*phi);};
    const dg::x::DVec vol = dg::create::volume( g3d);
    auto error = [&]( const dg::x::DVec& value, const dg::x::DVec& expected){
        dg::x::DVec diff( value);
        dg::blas1::axpby( 1., expected, -1., diff);
        return sqrt( dg::blas2::dot( diff, vol, diff));
    };
    SECTION( "One round")
    {
        const dg::x::DVec expected = dg::evaluate( [=]( double R, double Z, double phi){
            int idx = ((int)floor(phi/hz) - (int)p0 + Nz)%Nz;
            if( idx == 0)
                return binary(R,Z)*unary(0);
            return binary(R,Z)*( unary( idx*hz) + unary( -(int)(Nz-idx)*hz));
        }, g3d);
        CHECK( error( fa.evaluate( binary, unary, p0, 1), expected) < 1e-12);
        const dg::x::DVec free = dg::construct<dg::x::DVec>(
            dg::geo::fieldaligned_evaluate( g3d, bhat, binary, unary, p0, 1));
        CHECK( error( free, expected) < 1e-12);
    }
    SECTION( "No rounds")
    {
        const dg::x::DVec expected = dg::evaluate( [=]( double R, double Z, double phi){
            int idx = (int)floor(phi/hz) - (int)p0;
            return binary(R,Z)*unary( idx*hz);
        }, g3d);
        CHECK( error( fa.evaluate( binary, unary, p0, 0), expected) < 1e-12);
        const dg::x::DVec free = dg::construct<dg::x::DVec>(
            dg::geo::fieldaligned_evaluate( g3d, bhat, binary, unary, p0, 0));
        CHECK( error( free, expected) < 1e-12);
    }
}

#ifndef WITH_MPI // these functions do not exist in the MPI version
TEST_CASE( "Fieldaligned coarse grid functions")
{
    // In a toroidal field the plus and minus transformations just shift planes
    const unsigned n = 3, Nx = 16, Ny = 16;
    const dg::CylindricalGrid3d fine( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        n, Nx, Ny, 8, dg::NEU, dg::NEU, dg::PER);
    const dg::CylindricalGrid3d coarse( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        n, Nx, Ny, 4, dg::NEU, dg::NEU, dg::PER);
    const std::string method = GENERATE( "dg", "linear", "cubic", "linear-nearest",
        "cubic-nearest", "dg-nearest", "dg-linear", "linear-linear");
    INFO( "method "<<method);
    const dg::geo::TokamakMagneticField mag = dg::geo::createToroidalField( R_0);
    dg::geo::Fieldaligned<dg::aProductGeometry3d,dg::IDMatrix,dg::DVec> fa(
        dg::geo::createBHat(mag), fine, dg::NEU, dg::NEU, dg::geo::NoLimiter(),
        1e-8, 3, 3, -1, method, false);
    // Interpolation in the perpendicular plane is exact only for dg
    const double eps = method == "dg" ? 1e-13 : 1e-2;
    auto perp = [=]( double R, double Z){ return exp(R-R_0)*(1+Z*Z);};
    SECTION( "Interpolate is linear in phi")
    {
        // Plane i of the coarse grid is plane 2i of the fine grid (planes are
        // matched by index) and plane 2i+1 is the average of its neighbors
        const double hc = coarse.hz(), hf = fine.hz();
        auto modulate = [=]( double phi){ return 1. + cos(phi);};
        const dg::DVec in = dg::evaluate( [=]( double R, double Z, double phi){
            return perp(R,Z)*modulate(phi);}, coarse);
        dg::DVec out = fa.interpolate_from_coarse_grid( coarse, in);
        const dg::DVec expected = dg::evaluate( [=]( double R, double Z, double phi){
            unsigned k = floor( phi/hf);
            double left = (k/2 + 0.5)*hc;
            double alpha = (k%2)/2.;
            return perp(R,Z)*( (1.-alpha)*modulate(left) + alpha*modulate(left+hc));
        }, fine);
        const dg::DVec vol = dg::create::volume( fine);
        dg::blas1::axpby( 1., expected, -1., out);
        const double error = sqrt( dg::blas2::dot( out, vol, out)/dg::blas2::dot(
                expected, vol, expected));
        INFO( "Relative error "<<error);
        CHECK( error < eps);
    }
    SECTION( "Interpolate a constant")
    {
        const dg::DVec in = dg::evaluate( dg::one, coarse);
        dg::DVec out = fa.interpolate_from_coarse_grid( coarse, in);
        dg::blas1::axpby( 1., 1., -1., out);
        CHECK( dg::blas1::dot( out, out) < 1e-20);
    }
    SECTION( "Integrate between planes")
    {
        const dg::ClonePtr<dg::aGeometry2d> perp_ptr = coarse.perp_grid();
        const dg::DVec in = dg::evaluate( perp, *perp_ptr);
        dg::DVec out( in);
        fa.integrate_between_coarse_grid( coarse, in, out);
        const dg::DVec vol = dg::create::volume( *perp_ptr);
        dg::blas1::axpby( 1., in, -1., out);
        const double error = sqrt( dg::blas2::dot( out, vol, out)/dg::blas2::dot(
                in, vol, in));
        INFO( "Relative error "<<error);
        CHECK( error < eps);
        // out may alias in
        dg::DVec alias( in), result( in);
        fa.integrate_between_coarse_grid( coarse, in, result);
        fa.integrate_between_coarse_grid( coarse, alias, alias);
        CHECK( alias == result);
    }
}
#endif // WITH_MPI

TEST_CASE( "Fieldaligned construction")
{
    // This test is taken from a failed MAST simulation which failed on construction
#ifdef WITH_MPI
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {false,false,true});
#endif
    const dg::x::CylindricalGrid3d g3d( 96.2339, 752.624, -682.386, 682.386, 0, 2.*M_PI,
        3, 32, 67, 32, dg::NEU, dg::NEU, dg::PER
#ifdef WITH_MPI
    , comm
#endif
    );
    double R_0 = 469.921339;
    const dg::geo::TokamakMagneticField mag = dg::geo::createToroidalField( R_0);
    auto bhat = dg::geo::createBHat(mag);
    std::string method = "linear-nearest";
    dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>  dsFA(
            bhat, g3d, dg::NEU, dg::NEU, dg::geo::NoLimiter(), 1e-8, 12, 12,
            -1, method, false);
    CHECK(true);
}

#ifdef WITH_MPI
// The MPI Fieldaligned must give bitwise the same results as the shared
// memory Fieldaligned on the global grid on the same backend, for any
// number of ranks and any domain decomposition

// Compare the local MPI result with the expected local values: the
// differences may be at most max_ulp units in the last place of the largest
// expected value on this rank (0 means bitwise equal). A different order of
// summation changes a sum by a few ulp of its largest terms, which can be
// many ulp of a small result
static void check_ulp( std::string what, const dg::HVec& local,
    const dg::HVec& expected, unsigned max_ulp)
{
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    double max_value = 0;
    for( unsigned i=0; i<expected.size(); i++)
        max_value = std::max( max_value, fabs( expected[i]));
    const double ulp = std::nextafter( max_value, 2*max_value+1) - max_value;
    unsigned mismatches = 0;
    double worst = 0;
    for( unsigned i=0; i<local.size(); i++)
    {
        if( std::memcmp( &local[i], &expected[i], sizeof(double)) != 0)
            mismatches++;
        worst = std::max( worst, fabs( local[i] - expected[i])/ulp);
    }
    INFO( what<<" on rank "<<rank<<": "<<mismatches<<" of "<<local.size()
        <<" values differ, at most by "<<worst<<" ulp of the largest value "
        <<max_value<<" (allowed "<<max_ulp<<")");
    CHECK( worst <= max_ulp);
}


TEST_CASE( "MPI Fieldaligned agrees with shared memory Fieldaligned")
{
    int size;
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    // Distribute the ranks in the plane only, in z only, or in all directions
    std::array<int,3> dims = GENERATE( std::array<int,3>{0,0,1},
        std::array<int,3>{1,1,0}, std::array<int,3>{0,0,0});
    // n < 3 interpolates the fieldline end points with cubic polynomials
    unsigned n = GENERATE( 2, 3);
    std::string method = GENERATE( "dg", "linear-nearest", "dg-linear",
        "linear-equi");
    INFO( "Number of ranks "<<size<<", dims "<<dims[0]<<" "<<dims[1]<<" "
        <<dims[2]<<", n "<<n<<", method "<<method);
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD,
        {dims[0], dims[1], dims[2]}, {0,0,1});
    // The cell numbers are not multiples of the number of ranks, so that
    // the local grids have different sizes; 8 planes so that every process
    // has one with up to 8 processes in z; mx and my are multiples of n as
    // the projection methods other than dg require
    const unsigned Nx = 7, Ny = 11, Nz = 8, mx = 2*n, my = 2*n;
    dg::CylindricalMPIGrid3d g( R_0-a, R_0+a, -a, a, 0, 2.*M_PI, n, Nx, Ny, Nz,
        dg::NEU, dg::NEU, dg::PER, comm);
    dg::CylindricalGrid3d g_global( R_0-a, R_0+a, -a, a, 0, 2.*M_PI, n, Nx, Ny,
        Nz, dg::NEU, dg::NEU, dg::PER);
    auto bhat = dg::geo::createBHat( dg::geo::createCircularField( R_0, 20., a));

    dg::geo::Fieldaligned<dg::aProductMPIGeometry3d, dg::MIDMatrix, dg::MDVec>
        fa( bhat, g, dg::NEU, dg::NEU, dg::geo::NoLimiter(), 1e-8, mx, my,
            -1, method, false);
    dg::geo::Fieldaligned<dg::aProductGeometry3d, dg::IDMatrix, dg::DVec>
        fa_global( bhat, g_global, dg::NEU, dg::NEU, dg::geo::NoLimiter(),
            1e-8, mx, my, -1, method, false);

    auto local_part = [&]( const dg::DVec& global){
        return dg::global2local( dg::HVec( global), g).data();
    };
    // Fieldline integration
    check_ulp( "hbp", fa.hbp().data(), local_part( fa_global.hbp()), 0);
    check_ulp( "hbm", fa.hbm().data(), local_part( fa_global.hbm()), 0);
    check_ulp( "sqrtGp", fa.sqrtGp().data(), local_part( fa_global.sqrtGp()), 0);
    check_ulp( "sqrtGm", fa.sqrtGm().data(), local_part( fa_global.sqrtGm()), 0);
    check_ulp( "bphiP", fa.bphiP().data(), local_part( fa_global.bphiP()), 0);
    check_ulp( "bphiM", fa.bphiM().data(), local_part( fa_global.bphiM()), 0);
    // Fieldline integration in fieldaligned_evaluate
    auto binary = []( double R, double Z){ return sin( R)*cos( 3.*Z) + 2.;};
    auto unary = []( double s){ return exp( -s*s);};
    check_ulp( "fieldaligned_evaluate", dg::geo::fieldaligned_evaluate( g,
        bhat, binary, unary, 0, 2).data(), local_part( dg::DVec(
        dg::geo::fieldaligned_evaluate( g_global, bhat, binary, unary, 0, 2))),
        0);
    // evaluate applies the interpolation matrices
    check_ulp( "evaluate", fa.evaluate( binary, unary, 0, 2).data(),
        local_part( fa_global.evaluate( binary, unary, 0, 2)), 0);
    // Interpolation matrices
    dg::HVec global_f = dg::evaluate( []( double R, double Z, double phi){
        return sin( R)*cos( 3.*Z)*(1.5 + sin( phi)) + 1e-3*cos( 13.*R*Z);},
        g_global);
    dg::MDVec f = dg::global2local( global_f, g), result( f);
    dg::DVec gf = global_f, gresult( gf);
    for( auto which : { dg::geo::einsPlus, dg::geo::einsMinus,
        dg::geo::einsPlusT, dg::geo::einsMinusT, dg::geo::zeroPlus,
        dg::geo::zeroMinus})
    {
        INFO( "Matrix "<<which);
        fa( which, f, result);
        fa_global( which, gf, gresult);
        check_ulp( "matrix", result.data(), local_part( gresult), 0);
    }
}
#endif //WITH_MPI
