
#include <iostream>
#include "elliptic.h"
#include "catch2/catch_all.hpp"

double func1d ( double x) { return sin(2*M_PI*x);}
double lapM_func1d ( double x) { return 4*M_PI*M_PI*sin(2*M_PI*x);}

double func2d ( double x, double y) { return sin(2*M_PI*x)*sin(2*M_PI*y);}
double lapM_func2d ( double x, double y) { return 2*4*M_PI*M_PI*(sin(2*M_PI*x)*sin(2*M_PI*y));}

double func3d ( double x, double y, double z) { return sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);}
double lapM_func3d ( double x, double y, double z) { return 3*4*M_PI*M_PI*(sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z));}

TEST_CASE( "Elliptic1d")
{
    unsigned n = 19, N = 4;
    dg::Grid1d g1d( 0,1, n, N, dg::PER);
    const dg::DVec x = dg::evaluate( func1d, g1d);
    const dg::DVec w1d = dg::create::weights( g1d);
    const dg::DVec ysol = dg::evaluate( lapM_func1d, g1d);
    SECTION( "Construction")
    {
        const dg::Elliptic1d<dg::Grid1d, dg::DMatrix, dg::DVec> elliptic(g1d);
        CHECK( elliptic.get_jfactor() == 1);
        const auto& sigma = elliptic.get_chi();
        CHECK( dg::blas1::dot( sigma, sigma) == g1d.size());
    }
    SECTION( "Set and get")
    {
        dg::Elliptic1d<dg::Grid1d, dg::DMatrix, dg::DVec> elliptic(g1d);
        elliptic.set_jfactor( 10);
        CHECK( elliptic.get_jfactor() == 10);
        elliptic.set_chi( 4.0);
        const auto& sigma = elliptic.get_chi();
        CHECK( dg::blas1::dot( sigma, sigma) == 16*g1d.size());
    }

    SECTION( "Convergence")
    {
        auto dir = GENERATE( as<dg::direction>{}, dg::forward, dg::backward, dg::centered);
        const dg::Elliptic1d<dg::Grid1d, dg::DMatrix, dg::DVec> elliptic(g1d, dir);
        dg::DVec y(x);
        elliptic.symv( 1., x, 1., 2., 0, y);
        dg::blas1::axpby( 1., ysol, -1./2., y);
        CHECK ( dg::blas2::dot( y, w1d, y)  < 1e-16*dg::blas2::dot( x, w1d, x));
    }
    SECTION( "Convergence with sigma")
    {
        auto dir = GENERATE( as<dg::direction>{}, dg::forward, dg::backward, dg::centered);
        dg::Elliptic1d<dg::Grid1d, dg::DMatrix, dg::DVec> elliptic(g1d, dir);
        elliptic.set_chi( 10);
        dg::DVec y(x);
        elliptic.symv( 1., x, 0, y);
        dg::blas1::axpby( 1., ysol, -1./10., y);
        CHECK ( dg::blas2::dot( y, w1d, y)  < 1e-16*dg::blas2::dot( x, w1d, x));
    }
}
TEST_CASE( "Elliptic2d")
{
    unsigned nx = 19, ny = 17, Nx = 4, Ny = 5;
    dg::CartesianGrid2d g2d( dg::Grid1d{ 0, 1, nx, Nx}, dg::Grid1d{ 0, 1, ny, Ny} );
    const dg::DVec x = dg::evaluate( func2d, g2d);
    const dg::DVec w2d = dg::create::weights( g2d);
    const dg::DVec ysol = dg::evaluate( lapM_func2d, g2d);
    SECTION( "Construction")
    {
        const dg::Elliptic2d<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> elliptic(g2d);
        CHECK( elliptic.get_jfactor() == 1);
        const auto& sigma = elliptic.get_sigma();
        CHECK( dg::blas1::dot( sigma, sigma) == g2d.size());
    }
    SECTION( "Set and get")
    {
        dg::Elliptic2d<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> elliptic(g2d);
        elliptic.set_jfactor( 10);
        CHECK( elliptic.get_jfactor() == 10);
        elliptic.set_chi( 4.0);
        const auto& sigma = elliptic.get_sigma();
        CHECK( dg::blas1::dot( sigma, sigma) == 16*g2d.size());
    }

    SECTION( "Convergence")
    {
        auto dir = GENERATE( as<dg::direction>{}, dg::forward, dg::backward, dg::centered);
        const dg::Elliptic2d<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> elliptic(g2d, dir);
        dg::DVec y(x);
        elliptic.symv( 1., x, 1., 2., 0, y);
        dg::blas1::axpby( 1., ysol, -1./2., y);
        CHECK ( dg::blas2::dot( y, w2d, y)  < 1e-16*dg::blas2::dot( x, w2d, x));
    }
    SECTION( "Convergence with sigma")
    {
        auto dir = GENERATE( as<dg::direction>{}, dg::forward, dg::backward, dg::centered);
        dg::Elliptic2d<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> elliptic(g2d, dir);
        elliptic.set_chi( 10);
        dg::DVec y(x);
        elliptic.symv( 1., x, 0, y);
        dg::blas1::axpby( 1., ysol, -1./10., y);
        CHECK ( dg::blas2::dot( y, w2d, y)  < 1e-16*dg::blas2::dot( x, w2d, x));
    }
    SECTION( "Variation")
    {
        // sigma = lambda^2 |grad x|^2
        const dg::Elliptic2d<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> elliptic(g2d);
        const dg::DVec sol = dg::evaluate( []( double x, double y){
            double dx = 2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y);
            double dy = 2*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y);
            return 4.*(dx*dx + dy*dy);}, g2d);
        dg::DVec sigma(x);
        elliptic.variation( 2., x, sigma);
        dg::blas1::axpby( 1., sol, -1., sigma);
        CHECK( dg::blas2::dot( sigma, w2d, sigma) < 1e-16*dg::blas2::dot( sol, w2d, sol));
    }
    SECTION( "Jump weighting")
    {
        // For constant chi = c weighting the jumps with chi is the same as
        // the jump factor c; f jumps at the cell boundary x = 0.5
        const dg::DVec f = dg::evaluate( []( double x, double y){
            return (x < 0.5 ? 1. : 0.)*sin(2*M_PI*y);}, g2d);
        const double c = 3.;
        dg::Elliptic2d<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> weighted(g2d,
            dg::forward, 1., true), unweighted( g2d, dg::forward, c, false);
        CHECK( weighted.get_jump_weighting());
        CHECK( not unweighted.get_jump_weighting());
        weighted.set_chi( c);
        unweighted.set_chi( c);
        dg::DVec yw( f), yu( f);
        weighted.symv( f, yw);
        unweighted.symv( f, yu);
        dg::blas1::axpby( 1., yu, -1., yw);
        CHECK( dg::blas2::dot( yw, w2d, yw) < 1e-24*dg::blas2::dot( yu, w2d, yu));
        // and the jump term does not vanish
        unweighted.set_jfactor( 1.);
        unweighted.symv( f, yw);
        dg::blas1::axpby( 1., yu, -1., yw);
        CHECK( dg::blas2::dot( yw, w2d, yw) > 1e-8*dg::blas2::dot( yu, w2d, yu));
    }
}
TEST_CASE( "Elliptic3d")
{
    unsigned nx = 19, ny = 17, Nx = 4, Ny = 5, nz = 18, Nz = 6;
    dg::CartesianGrid3d g3d( dg::Grid1d{ 0, 1, nx, Nx}, dg::Grid1d{ 0, 1, ny, Ny}, dg::Grid1d{0,1,nz,Nz} );
    const dg::DVec x = dg::evaluate( func3d, g3d);
    const dg::DVec w3d = dg::create::weights( g3d);
    const dg::DVec ysol = dg::evaluate( lapM_func3d, g3d);
    SECTION( "Construction")
    {
        const dg::Elliptic3d<dg::CartesianGrid3d, dg::DMatrix, dg::DVec> elliptic(g3d);
        CHECK( elliptic.get_jfactor() == 1);
        const auto& sigma = elliptic.get_sigma();
        CHECK( dg::blas1::dot( sigma, sigma) == g3d.size());
    }
    SECTION( "Set and get")
    {
        dg::Elliptic3d<dg::CartesianGrid3d, dg::DMatrix, dg::DVec> elliptic(g3d);
        elliptic.set_jfactor( 10);
        CHECK( elliptic.get_jfactor() == 10);
        elliptic.set_chi( 4.0);
        const auto& sigma = elliptic.get_sigma();
        CHECK( dg::blas1::dot( sigma, sigma) == 16*g3d.size());
    }

    SECTION( "Convergence")
    {
        auto dir = GENERATE( as<dg::direction>{}, dg::forward, dg::backward, dg::centered);
        const dg::Elliptic3d<dg::CartesianGrid3d, dg::DMatrix, dg::DVec> elliptic(g3d, dir);
        dg::DVec y(x);
        elliptic.symv( 1., x, 1., 2., 0, y);
        dg::blas1::axpby( 1., ysol, -1./2., y);
        CHECK ( dg::blas2::dot( y, w3d, y)  < 1e-16*dg::blas2::dot( x, w3d, x));
    }
    SECTION( "Convergence with sigma")
    {
        auto dir = GENERATE( as<dg::direction>{}, dg::forward, dg::backward, dg::centered);
        dg::Elliptic3d<dg::CartesianGrid3d, dg::DMatrix, dg::DVec> elliptic(g3d, dir);
        elliptic.set_chi( 10);
        dg::DVec y(x);
        elliptic.symv( 1., x, 0, y);
        dg::blas1::axpby( 1., ysol, -1./10., y);
        CHECK ( dg::blas2::dot( y, w3d, y)  < 1e-16*dg::blas2::dot( x, w3d, x));
    }
    SECTION( "Compute in 2d")
    {
        // In 2d mode Elliptic3d is identical to Elliptic2d (also for nz > 1)
        // f jumps at the cell boundary z = 0.5 such that the jump terms in z
        // do not vanish
        const dg::DVec f = dg::evaluate( []( double x, double y, double z){
            return sin(2*M_PI*x)*sin(2*M_PI*y)*( z < 0.5 ? 1. : 0.);}, g3d);
        bool jump_weighting = GENERATE( false, true);
        INFO( "Jump weighting "<<jump_weighting);
        dg::Elliptic3d<dg::CartesianGrid3d, dg::DMatrix, dg::DVec> elliptic3d(
            g3d, dg::forward, 1., jump_weighting);
        dg::Elliptic2d<dg::CartesianGrid3d, dg::DMatrix, dg::DVec> elliptic2d(
            g3d, dg::forward, 1., jump_weighting);
        elliptic3d.set_compute_in_2d( true);
        dg::DVec y3d(f), y2d(f);
        elliptic3d.symv( f, y3d);
        elliptic2d.symv( f, y2d);
        dg::blas1::axpby( 1., y2d, -1., y3d);
        CHECK( dg::blas2::dot( y3d, w3d, y3d) < 1e-28*dg::blas2::dot( y2d, w3d, y2d));
    }
}

TEST_CASE( "Boundary conditions and directions")
{
    // Elliptic converts boundary conditions and directions with enums.h
    for( auto bc : { dg::PER, dg::DIR, dg::DIR_NEU, dg::NEU_DIR, dg::NEU})
    {
        CHECK( dg::str2bc( dg::bc2str( bc)) == bc);
        CHECK( dg::inverse( dg::inverse( bc)) == bc);
    }
    for( auto dir : { dg::forward, dg::backward, dg::centered})
    {
        CHECK( dg::str2direction( dg::direction2str( dir)) == dir);
        CHECK( dg::inverse( dg::inverse( dir)) == dir);
    }
}
