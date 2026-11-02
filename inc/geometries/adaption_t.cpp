#include <iostream>

#include "dg/algorithm.h"
#include "adaption.h"
#include "guenter.h"
#include "catch2/catch_all.hpp"

// 4th order central differences
static double diffX( const dg::geo::CylindricalFunctor& f, double x, double y)
{
    const double h = 1e-3;
    return (-f(x+2*h,y)+8*f(x+h,y)-8*f(x-h,y)+f(x-2*h,y))/(12*h);
}
static double diffY( const dg::geo::CylindricalFunctor& f, double x, double y)
{
    const double h = 1e-3;
    return (-f(x,y+2*h)+8*f(x,y+h)-8*f(x,y-h)+f(x,y-2*h))/(12*h);
}

TEST_CASE( "Adaption")
{
    const dg::geo::CylindricalFunctorsLvl2 psi = dg::geo::guenter::createPsip( 10);
    const double x = 10.3, y = 0.4, eps = 1e-6;
    SECTION( "Inverse gradient")
    {
        auto inv = dg::geo::make_NablaPsiInvCollective( psi);
        CHECK_THAT( inv.f()(x,y), Catch::Matchers::WithinRel( 1./hypot(
                psi.dfx()(x,y), psi.dfy()(x,y)), eps));
        CHECK_THAT( inv.dfx()(x,y), Catch::Matchers::WithinRel( diffX( inv.f(), x, y), eps));
        CHECK_THAT( inv.dfy()(x,y), Catch::Matchers::WithinRel( diffY( inv.f(), x, y), eps));
    }
    SECTION( "Liseikin monitor metric")
    {
        auto chi = dg::geo::make_LiseikinCollective( psi, 2.3, 0.7);
        CHECK_THAT( chi.xx()(x,y)*chi.yy()(x,y) - chi.xy()(x,y)*chi.xy()(x,y),
                Catch::Matchers::WithinRel( 1., eps));
        CHECK_THAT( chi.divX()(x,y), Catch::Matchers::WithinRel(
                diffX( chi.xx(), x, y) + diffY( chi.xy(), x, y), eps));
        CHECK_THAT( chi.divY()(x,y), Catch::Matchers::WithinRel(
                diffX( chi.xy(), x, y) + diffY( chi.yy(), x, y), eps));
    }
}
