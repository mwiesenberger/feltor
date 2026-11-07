#include <iostream>

#include "polarization.h"

#include "catch2/catch_all.hpp"

using Container = dg::DVec;
using Geometry = dg::CartesianGrid2d;

TEST_CASE( "PolCharge")
{
    // For f = sin(x) sin(y) on a periodic grid -Lap f = 2 f and the
    // Helmholtz operator 1 + alpha Lap has the Eigenvalue 1 - 2 alpha
    Geometry grid( 0., 2.*M_PI, 0, 2.*M_PI, 3, 16, 16, dg::PER, dg::PER);
    const Container w2d = dg::create::weights( grid);
    const Container f = dg::evaluate( []( double x, double y){
            return sin(x)*sin(y);}, grid);
    Container y( f);
    const double alpha = -0.5;
    const std::vector<double> eps = {1e-10, 1e-10, 1e-10};
    auto ratio = [&]( const Container& x){
        return dg::blas2::dot( f, w2d, x)/dg::blas2::dot( f, w2d, f);
    };
    SECTION( "Copy of ff")
    {
        // The copy must not refer to the operators of the original
        dg::mat::PolCharge<Geometry, dg::DMatrix, Container> copy;
        {
            dg::mat::PolCharge<Geometry, dg::DMatrix, Container> pol( alpha,
                    eps, grid, dg::centered, 1., "ff");
            copy = pol;
        }
        dg::blas2::symv( copy, f, y);
        INFO( "<f, P f>/<f,f> "<<ratio(y));
        // up to the discretization error of the Laplacian
        CHECK( fabs( ratio(y) - 2./(1.-2.*alpha)) < 1e-5);
    }
}
