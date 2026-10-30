#include <iostream>
#include <functional>

#include "backend/typedefs.h"
#include "elliptic.h"
#include "multigrid.h"

#include "catch2/catch_all.hpp"

static double sol( double x, double y) { return sin(x)*sin(y);}
static double rhs( double x, double y) { return 2.*sin(x)*sin(y);}

TEST_CASE( "Multigrid")
{
    // All multigrid methods must reproduce the PCG solution on the fine grid
    const unsigned n = 3, Nx = 32, Ny = 32;
    const double eps = 1e-12;
    const dg::CartesianGrid2d grid( 0, 2.*M_PI, 0, 2.*M_PI, n, Nx, Ny, dg::DIR, dg::DIR);
    const dg::DVec w2d = dg::create::weights( grid);
    const dg::DVec b = dg::evaluate( rhs, grid);
    dg::Elliptic<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> pol( grid);
    dg::DVec x_ref( b.size(), 0.);
    dg::PCG<dg::DVec> pcg( x_ref, grid.size());
    pcg.solve( pol, x_ref, b, pol.precond(), pol.weights(), eps);
    const double nrm_ref = dg::blas2::dot( x_ref, w2d, x_ref);
    // the reference solution is the discrete solution
    dg::DVec error = dg::evaluate( sol, grid);
    dg::blas1::axpby( 1., x_ref, -1., error);
    REQUIRE( dg::blas2::dot( error, w2d, error) < 1e-8*nrm_ref);

    unsigned stages = GENERATE( 1u, 3u);
    INFO( "Number of stages "<<stages);
    dg::NestedGrids<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> nested( grid,
        stages);
    std::vector<dg::Elliptic<dg::CartesianGrid2d, dg::DMatrix, dg::DVec>> ops(
        stages);
    for( unsigned u=0; u<stages; u++)
        ops[u].construct( nested.grid(u));
    dg::DVec x( b.size(), 0.);
    SECTION( "NestedGrids")
    {
        REQUIRE( nested.stages() == stages);
        // Each stage halves the number of cells
        for( unsigned u=0; u<stages; u++)
        {
            CHECK( nested.grid(u).Nx() == Nx >> u);
            CHECK( nested.grid(u).Ny() == Ny >> u);
        }
        // The projection preserves constants
        std::vector<dg::DVec> ones = nested.project( dg::DVec( b.size(), 1.));
        REQUIRE( ones.size() == stages);
        for( unsigned u=0; u<stages; u++)
        {
            dg::DVec diff = ones[u];
            dg::blas1::plus( diff, -1.);
            CHECK( diff.size() == nested.grid(u).size());
            CHECK( dg::blas1::dot( diff, diff) < 1e-28*diff.size());
        }
    }
    SECTION( "MultigridCG2d")
    {
        dg::MultigridCG2d<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> multigrid(
            grid, stages);
        multigrid.set_benchmark( false);
        std::vector<unsigned> number = multigrid.solve( ops, x, b, eps);
        CHECK( number.size() == stages);
        dg::blas1::axpby( 1., x_ref, -1., x);
        CHECK( dg::blas2::dot( x, w2d, x) < 1e-12*nrm_ref);
    }
    SECTION( "fmg_solve")
    {
        // Use PCG as inverse operators on all stages
        std::vector<dg::PCG<dg::DVec>> pcgs( stages);
        std::vector<std::function<void( const dg::DVec&, dg::DVec&)>> inverse(
            stages);
        for( unsigned u=0; u<stages; u++)
        {
            pcgs[u].construct( nested.x(u), nested.grid(u).size());
            inverse[u] = [eps, &pcg = pcgs[u], &op = ops[u]]( const auto& y,
                auto& z)
            {
                pcg.solve( op, z, y, op.precond(), op.weights(), eps);
            };
        }
        std::vector<std::function<void( const dg::DVec&, dg::DVec&)>> down(
            inverse.begin(), inverse.end()-1);
        dg::fmg_solve( ops, x, b, down, inverse, nested, w2d, eps, 1);
        dg::blas1::axpby( 1., x_ref, -1., x);
        CHECK( dg::blas2::dot( x, w2d, x) < 1e-12*nrm_ref);
    }
}
