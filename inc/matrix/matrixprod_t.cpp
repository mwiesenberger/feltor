#include <iostream>

#include "matrixprod.h"

#include "catch2/catch_all.hpp"

using Container = dg::DVec;

TEST_CASE( "ProductMatrixFunction")
{
    // For a diagonal matrix A and a vector d both f(d, A) b and f(A, d) b are
    // the pointwise product f( d, a) b
    dg::CartesianGrid2d grid( 0., 2.*M_PI, 0, 2.*M_PI, 3, 16, 16, dg::PER, dg::PER);
    const Container w2d = dg::create::weights( grid);
    const Container a = dg::evaluate( []( double x, double y){
            return 2. + sin(x)*cos(y);}, grid);
    const Container d = dg::evaluate( []( double x, double y){
            return 1. + 0.5*cos(x+y);}, grid);
    const Container b = dg::evaluate( []( double x, double y){
            return 1.5 + cos(x+0.5)*sin(2.*y);}, grid);
    auto diagonal = [&a]( const auto& y, auto& z){ dg::blas1::pointwiseDot( a, y, z);};
    auto op = []DG_DEVICE( double dd, double lambda){ return exp( -dd*lambda);};
    Container x(b), solution(b), error(b);
    dg::blas1::evaluate( solution, dg::equals(), []DG_DEVICE( double dd,
            double aa, double bb){ return exp( -dd*aa)*bb;}, d, a, b);
    auto rel_error = [&]()
    {
        dg::blas1::axpby( 1., solution, -1., x, error);
        return sqrt( dg::blas2::dot( w2d, error)/dg::blas2::dot( w2d, solution));
    };
    dg::mat::ProductMatrixFunction<Container> prod( w2d, 200);
    prod.set_benchmark( false);
    SECTION( "apply")
    {
        unsigned iter = prod.apply( x, op, d, diagonal, b, w2d, 1e-10);
        INFO( "Iterations "<<iter<<" error "<<rel_error());
        CHECK( iter > 1);
        CHECK( rel_error() < 1e-9);
    }
    SECTION( "compute_vlcl with a given tridiagonalisation")
    {
        auto func = dg::mat::make_FuncEigen_Te1( [&]( double t){ return op( 1., t);});
        auto T = prod.lanczos().tridiag( func, diagonal, b, w2d, 1e-10, 1.,
                "universal", 1., 2);
        prod.compute_vlcl( op, d, diagonal, T, x, b, prod.lanczos().get_bnorm());
        CHECK( rel_error() < 1e-9);
    }
    SECTION( "Tridiagonal matrix of size 1")
    {
        // b is an Eigenvector of a constant diagonal, Lanczos stops after 1 iteration
        const Container two( b.size(), 2.);
        auto constant = [&two]( const auto& y, auto& z){ dg::blas1::pointwiseDot( two, y, z);};
        unsigned iter = prod.apply( x, op, d, constant, b, w2d, 1e-10);
        CHECK( iter == 1);
        dg::blas1::evaluate( solution, dg::equals(), []DG_DEVICE( double dd,
                double bb){ return exp( -2.*dd)*bb;}, d, b);
        CHECK( rel_error() < 1e-14);
    }
    // apply_adjoint is not tested: it is documented not to converge
}
