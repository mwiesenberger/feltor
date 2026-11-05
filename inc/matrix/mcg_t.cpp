#include <iostream>

#include "lanczos.h"
#include "mcg.h"

#include "catch2/catch_all.hpp"

using Container = dg::DVec;

TEST_CASE( "MCGFuncEigen")
{
    // Multiplication with d(x,y) has the analytic function f(d) b
    dg::CartesianGrid2d grid( 0., 2.*M_PI, 0, 2.*M_PI, 3, 16, 16, dg::PER, dg::PER);
    const Container w2d = dg::create::weights( grid);
    const Container d = dg::evaluate( []( double x, double y){
            return 2. + sin(x)*cos(y);}, grid);
    const Container b = dg::evaluate( []( double x, double y){
            return 1.5 + cos(x+0.5)*sin(2.*y);}, grid);
    auto diagonal = [&d]( const auto& y, auto& z){ dg::blas1::pointwiseDot( d, y, z);};
    auto func = []( double t){ return exp(-t);};
    Container solution(b), x(b), error(b);
    dg::blas1::evaluate( solution, dg::equals(), [&func]( double dd, double vv){
            return func(dd)*vv;}, d, b);
    auto rel_error = [&]( const Container& sol, const Container& y)
    {
        dg::blas1::axpby( 1., sol, -1., y, error);
        return sqrt( dg::blas2::dot( w2d, error)/dg::blas2::dot( w2d, sol));
    };

    dg::mat::MCGFuncEigen<Container> mcg( w2d, 500);
    mcg( x, func, diagonal, b, w2d, 1e-12, 1., 1.);
    INFO( "First call");
    CHECK( rel_error( solution, x) < 1e-12);

    SECTION( "Same as Lanczos")
    {
        // Equivalent to Lanczos with the residual error norm
        dg::mat::UniversalLanczos<Container> lanczos( w2d, 500);
        Container y(b);
        lanczos.solve( y, dg::mat::make_FuncEigen_Te1( func), diagonal, b,
                w2d, 1e-12, 1., "residual");
        CHECK( rel_error( y, x) < 1e-12);
    }
    SECTION( "Later calls")
    {
        // Results must not depend on previous calls
        mcg( x, func, diagonal, b, w2d, 1e-12, 1., 1.);
        INFO( "Second call");
        CHECK( rel_error( solution, x) < 1e-12);
        mcg( x, func, diagonal, b, w2d, 1e-4, 1., 1.);
        INFO( "Fewer iterations");
        CHECK( rel_error( solution, x) < 1e-4);
        mcg( x, func, diagonal, b, w2d, 1e-12, 1., 1.);
        INFO( "More iterations than the previous call");
        CHECK( rel_error( solution, x) < 1e-12);
    }
}
