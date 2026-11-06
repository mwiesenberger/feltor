#include <iostream>

#include "matrixsqrt.h"
#include "outer.h"

#include "catch2/catch_all.hpp"

using Container = dg::DVec;

TEST_CASE( "LaplaceDecomposition")
{
    // f(-Lap) b from the Eigen-decomposition of the 1d Laplacians equals the
    // Lanczos approximation of f applied to dg::Elliptic
    auto Ns = GENERATE( std::array<unsigned,2>{12,6}, std::array<unsigned,2>{6,12});
    INFO( "Nx "<<Ns[0]<<" Ny "<<Ns[1]);
    dg::CartesianGrid2d grid( 0., 2.*M_PI, 0, 2.*M_PI, 3, Ns[0], Ns[1], dg::PER, dg::PER);
    const Container w2d = dg::create::weights( grid);
    const Container b = dg::evaluate( []( double x, double y){
            return 1.5 + cos(x+0.5)*sin(2.*y);}, grid);
    Container x( b), y( b), error( b);

    dg::Grid2d g2d( 0., 2.*M_PI, 0, 2.*M_PI, 3, Ns[0], Ns[1], dg::PER, dg::PER);
    dg::mat::LaplaceDecomposition<Container> decomposition( g2d, dg::PER, dg::PER, dg::centered);
    decomposition.matrix_function( x, []( double t){ return exp(-t);}, b, 1e-12);

    dg::Elliptic<dg::CartesianGrid2d, dg::DMatrix, Container> ell( {grid, dg::centered});
    dg::mat::MatrixFunction<Container> func( ell, w2d, 1e-12);
    func.set_benchmark( false);
    func( []( double t){ return exp(-t);}, b, y);

    dg::blas1::axpby( 1., y, -1., x, error);
    double rel = sqrt( dg::blas2::dot( w2d, error)/dg::blas2::dot( w2d, y));
    INFO( "Relative difference "<<rel);
    CHECK( rel < 1e-10);
}
