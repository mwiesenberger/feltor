
#include <iostream>
#include "functors.h"

#include "lanczos.h"

#include "catch2/catch_all.hpp"

const double lx = 2.*M_PI;
const double ly = 2.*M_PI;
const dg::bc bcx = dg::DIR;
const dg::bc bcy = dg::DIR;

using Matrix = dg::DMatrix;
using Container = dg::DVec;

TEST_CASE( "Lanczos")
{
    unsigned n = 3, Nx = 16, Ny = 16;
    dg::CartesianGrid2d grid( 0., lx, 0, ly, n, Nx, Ny, bcx, bcy);
    const Container w2d = dg::create::weights( grid);
    const Container rnd = dg::construct<Container>(
            dg::create::uniform_real_distribution( 0, 1, grid));

    dg::mat::UniversalLanczos<Container> lanczos( w2d, 2000);
    SECTION( "Dirichlet")
    {
        grid.set_bcs( {dg::DIR, dg::DIR});
        dg::Elliptic<dg::CartesianGrid2d, Matrix, Container> ell( {grid, dg::centered, 1.0});
        //lanczos.set_verbose(true);
        auto T = lanczos.tridiag( ell, rnd, w2d);
        auto extremeEVs = dg::mat::compute_extreme_EV( T);
        // Eigenvalues computed with pyfeltor and scipy
        CHECK( extremeEVs[0] - 0.5000000674740441 < 1e-3);
        CHECK( 895.3216230750616 - extremeEVs[1] < 1e-3*895);
    }
    SECTION( "NEU")
    {
        grid.set_bcs( {dg::NEU, dg::NEU});
        dg::Elliptic<dg::CartesianGrid2d, Matrix, Container> ell( {grid, dg::centered, 1.0});
        //lanczos.set_verbose(true);
        auto T = lanczos.tridiag( ell, rnd, w2d);
        auto extremeEVs = dg::mat::compute_extreme_EV( T);

        CHECK( extremeEVs[0] - 0.0 < 1e-10*897);
        CHECK( 897.6777549436483 - extremeEVs[1] < 1e-3*897);
    }
    SECTION( "PER")
    {
        grid.set_bcs( {dg::PER, dg::PER});
        dg::Elliptic<dg::CartesianGrid2d, Matrix, Container> ell( {grid, dg::centered, 1.0});
        //lanczos.set_verbose(true);
        auto T = lanczos.tridiag( ell, rnd, w2d );
        auto extremeEVs = dg::mat::compute_extreme_EV( T);

        CHECK( extremeEVs[0] - 0.0 < 1e-10*897);
        CHECK( 896.6482950890147 - extremeEVs[1] < 1e-3*897);

        auto Tmax = lanczos.tridiag( ell, rnd, w2d, 1e-4, 1, "compute_max_EV" );
        extremeEVs = dg::mat::compute_extreme_EV( Tmax);

        CHECK( 896.6482950890147 - extremeEVs[1] < 1e-3*897);
        // the max EV needs fewer iterations!
        CHECK( Tmax.size() < T.size());
    }
}

TEST_CASE( "Lanczos solve")
{
    // Whenever the iteration stops, x must be computed from the final T
    unsigned n = 3, Nx = 16, Ny = 16;
    dg::CartesianGrid2d grid( 0., lx, 0, ly, n, Nx, Ny, dg::PER, dg::PER);
    const Container w2d = dg::create::weights( grid);
    const Container one = dg::evaluate( dg::one, grid);
    Container x(one), error(one);

    auto func = []( double t){ return exp(-t);};
    auto expT = dg::mat::make_FuncEigen_Te1( func);
    auto rel_error = [&]( const Container& solution)
    {
        dg::blas1::axpby( 1., solution, -1., x, error);
        return sqrt( dg::blas2::dot( w2d, error)/dg::blas2::dot( w2d, solution));
    };
    dg::Elliptic<dg::CartesianGrid2d, Matrix, Container> ell( {grid, dg::centered, 1.0});
    // Multiplication with d(x,y) has an analytic function
    auto make_diagonal = [&]( const Container& d)
    {
        return [&d]( const auto& y, auto& z){ dg::blas1::pointwiseDot( d, y, z);};
    };
    const Container smooth = dg::evaluate( []( double x, double y){
            return 2. + sin(x)*cos(y);}, grid);
    const Container b = dg::evaluate( []( double x, double y){
            return 1.5 + cos(x+0.5)*sin(2.*y);}, grid);
    // exp(-d) v
    Container solution( b);
    dg::mat::UniversalLanczos<Container> lanczos( w2d, 500);
    SECTION( "Breakdown in the first iteration, first call")
    {
        // A one = 0, there is no T for which the universal error was computed
        unsigned iter = lanczos.solve( x, expT, ell, one, w2d, 1e-8);
        CHECK( iter == 1);
        CHECK( rel_error( one) < 1e-10);
    }
    SECTION( "Breakdown in the first iteration, later call")
    {
        auto diagonal = make_diagonal( smooth);
        lanczos.solve( x, expT, diagonal, b, w2d, 1e-8);
        unsigned iter = lanczos.solve( x, expT, ell, one, w2d, 1e-8);
        CHECK( iter == 1);
        CHECK( rel_error( one) < 1e-10);
    }
    SECTION( "Breakdown between two checks")
    {
        // Three distinct Eigenvalues: T has size 3
        const Container steps = dg::evaluate( []( double x, double){
                return x < 2. ? 1. : ( x < 4. ? 2. : 3.);}, grid);
        auto diagonal = make_diagonal( steps);
        unsigned iter = lanczos.solve( x, expT, diagonal, b, w2d, 1e-12,
                1., "universal", 1., 1);
        CHECK( iter == 3);
        dg::blas1::evaluate( solution, dg::equals(), []DG_DEVICE( double dd,
                double vv){ return exp(-dd)*vv;}, steps, b);
        CHECK( rel_error( solution) < 1e-10);
    }
    SECTION( "Size one T with a function that uses T.P")
    {
        auto linear = dg::mat::make_Linear_Te1( +1);
        lanczos.solve( x, linear, ell, one, w2d, 1e-8);
        CHECK( sqrt( dg::blas2::dot( w2d, x)) < 1e-10); // A one = 0
    }
    SECTION( "Eigenvalue error norm")
    {
        auto diagonal = make_diagonal( smooth);
        lanczos.solve( x, expT, diagonal, b, w2d, 1e-8, 1., "compute_extreme_EV");
        dg::blas1::evaluate( solution, dg::equals(), []DG_DEVICE( double dd,
                double vv){ return exp(-dd)*vv;}, smooth, b);
        CHECK( rel_error( solution) < 1e-8);
    }
    SECTION( "Universal error with q = 2 is symmetric")
    {
        // The ODE solve of sqrt(T) needs the extended T to be symmetric
        const Container wide = dg::evaluate( []( double x, double){
                return 1. + 4000.*x*x/lx/lx;}, grid);
        auto diagonal = make_diagonal( wide);
        unsigned number = 0;
        auto sqrtT = dg::mat::make_SqrtODE_Te1( +1,
                std::string("Dormand-Prince-7-4-5"), 1e-8, 1e-10, number);
        lanczos.solve( x, sqrtT, diagonal, b, w2d, 1e-6, 1., "universal", 1., 2);
        dg::blas1::evaluate( solution, dg::equals(), []DG_DEVICE( double dd,
                double vv){ return sqrt(dd)*vv;}, wide, b);
        CHECK( rel_error( solution) < 1e-5);
    }
    SECTION( "Function that vanishes on the spectrum")
    {
        // exp(-1000 t) underflows to 0 for all Eigenvalues in [1,3]
        auto zero = dg::mat::make_FuncEigen_Te1( []( double t){ return exp(-1000.*t);});
        auto diagonal = make_diagonal( smooth);
        lanczos.solve( x, zero, diagonal, b, w2d, 1e-8);
        CHECK( sqrt( dg::blas2::dot( w2d, x)) == 0);
    }
    SECTION( "Large right hand side")
    {
        // Must not stop before the first check of the universal error
        Container large(b);
        dg::blas1::scal( large, 1e20);
        auto diagonal = make_diagonal( smooth);
        unsigned iter = lanczos.solve( x, expT, diagonal, large, w2d, 1e-8);
        CHECK( iter > 1);
        dg::blas1::evaluate( solution, dg::equals(), []DG_DEVICE( double dd,
                double vv){ return exp(-dd)*vv;}, smooth, large);
        CHECK( rel_error( solution) < 1e-7);
    }
}
