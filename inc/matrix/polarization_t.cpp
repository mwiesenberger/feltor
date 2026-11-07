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
    SECTION( "Modes")
    {
        // "df": -Lap (1+alpha Lap)^{-1}, "ff": (1+alpha Lap)^{-1/2} (-Lap)
        // (1+alpha Lap)^{-1/2} (chi = 1), "ffO4": (1+alpha Lap)^{-1}
        // (-Lap - Lap Lap + 2 div div grad grad)(1+alpha Lap)^{-1} (chi = iota = 1)
        auto mode = GENERATE( as<std::string>{}, "df", "ff", "ffO4");
        auto a = GENERATE( 0., -0.5);
        INFO( "Mode "<<mode<<" alpha "<<a);
        dg::mat::PolCharge<Geometry, dg::DMatrix, Container> pol( a, eps, grid,
                dg::centered, 1., mode);
        dg::blas2::symv( pol, f, y);
        double expected = mode == "ffO4" ? 6./(1.-2.*a)/(1.-2.*a) : 2./(1.-2.*a);
        INFO( "<f, P f>/<f,f> "<<ratio(y)<<" expected "<<expected);
        // up to the discretization error of the Laplacian
        CHECK( fabs( ratio(y) - expected) < 1e-4*expected);
    }
    SECTION( "df with commute")
    {
        // (1+alpha Lap)^{-1} (-Lap) has the same Eigenvalue
        dg::mat::PolCharge<Geometry, dg::DMatrix, Container> pol( alpha, eps,
                grid, dg::centered, 1., "df", true);
        dg::blas2::symv( pol, f, y);
        CHECK( fabs( ratio(y) - 2./(1.-2.*alpha)) < 1e-5);
    }
    SECTION( "Not implemented")
    {
        CHECK_THROWS_AS( (dg::mat::PolCharge<Geometry, dg::DMatrix, Container>(
                alpha, eps, grid, dg::centered, 1., "typo")), dg::Error);
        dg::mat::PolCharge<Geometry, dg::DMatrix, Container> pol( alpha, eps,
                grid, dg::centered, 1., "ff", true);
        CHECK_THROWS_AS( dg::blas2::symv( pol, f, y), dg::Error);
    }
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
