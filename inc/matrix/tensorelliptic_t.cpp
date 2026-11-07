#include <iostream>

#include "tensorelliptic.h"

#include "catch2/catch_all.hpp"

using Container = dg::DVec;

TEST_CASE( "TensorElliptic")
{
    // Constant chi and iota on a periodic grid: for f = sin(x) sin(y)
    // (-div chi grad - Lap iota Lap + div div 2 iota grad grad) f = (2 chi + 4 iota) f
    dg::CartesianGrid2d grid( 0., 2.*M_PI, 0, 2.*M_PI, 3, 32, 32, dg::PER, dg::PER);
    const Container w2d = dg::create::weights( grid);
    const Container f = dg::evaluate( []( double x, double y){
            return sin(x)*sin(y);}, grid);
    Container y( f);
    auto rel_error = [&]( const Container& x, const Container& solution)
    {
        Container error( x);
        dg::blas1::axpby( 1., solution, -1., x, error);
        return sqrt( dg::blas2::dot( w2d, error)/dg::blas2::dot( w2d, solution));
    };
    dg::mat::TensorElliptic<dg::CartesianGrid2d, dg::DMatrix, Container> tensor(
            grid, dg::centered, 1.);
    const double chi = 0.7, iota = 0.3;
    tensor.set_chi( Container( f.size(), chi));
    tensor.set_iota( Container( f.size(), iota));
    SECTION( "symv")
    {
        // Applying a dG second derivative to a second derivative is accurate
        // only in the weak sense, so test with a smooth function h
        dg::blas2::symv( tensor, f, y);
        const Container h = dg::evaluate( []( double x, double y){
                return sin(x)*sin(y) + cos(x)*sin(y);}, grid);
        double weak = dg::blas2::dot( h, w2d, y)/dg::blas2::dot( h, w2d, f);
        INFO( "<h, T f> / <h, f> "<<weak<<" expected "<<2.*chi+4.*iota);
        CHECK( fabs( weak - (2.*chi + 4.*iota)) < 1e-4);
    }
    SECTION( "variation")
    {
        // -chi/2 ( |grad phi|^2 + alpha chi ( |grad grad phi|^2 - (Lap phi)^2/2))
        // with |grad grad phi|^2 - (Lap phi)^2/2 = 2 cos^2 x cos^2 y
        const double alpha = 0.4, chiv = 1.3;
        const Container chivec( f.size(), chiv);
        tensor.variation( f, alpha, chivec, y);
        const Container solution = dg::evaluate( [&]( double x, double y){
                double grad2 = cos(x)*cos(x)*sin(y)*sin(y) + sin(x)*sin(x)*cos(y)*cos(y);
                return -chiv/2.*( grad2 + alpha*chiv*2.*cos(x)*cos(x)*cos(y)*cos(y));
            }, grid);
        INFO( "Relative error "<<rel_error( y, solution));
        CHECK( rel_error( y, solution) < 1e-3);
    }
}
