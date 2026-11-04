#include <iostream>
#include <iomanip>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "dg/algorithm.h"
#include "magnetic_field.h"
#include "testfunctors.h"
#include "ds.h"
#include "toroidal.h"
#include "catch2/catch_all.hpp"

const double R_0 = 10;
const double I_0 = 20; //q factor at r=1 is I_0/R_0
const double a  = 1; //small radius

TEST_CASE( "DS")
{
#ifdef WITH_MPI
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {false,false,true});
#endif
    INFO( "# Test the parallel derivative DS in cylindrical coordinates for \
circular flux surfaces with DIR and NEU boundary conditions.");
    std::string method = GENERATE( "dg", "linear", "cubic");
    unsigned mm = GENERATE( 1, 10);
    dg::bc bcx = GENERATE( dg::NEU, dg::DIR);
    unsigned n = 3, Nx = 26, Ny = 26, Nz = 20, mx[2] = {mm,mm}, max_iter = 1e4;

    INFO( "Combination:\n"
              <<"n:  "<<n<<"\n"
              <<"Nx: "<<Nx<<"\n"
              <<"Ny: "<<Ny<<"\n"
              <<"Nz: "<<Nz<<"\n"
              <<"mx: "<<mx[0]<<"\n"
              <<"my: "<<mx[1]<<"\n"
              <<"method: "<< method << "\n"
              <<"bc: "<<dg::bc2str(bcx)<<"\n"
              );
    const dg::x::CylindricalGrid3d g3d( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        n, Nx, Ny, Nz, bcx, bcx, dg::PER
#ifdef WITH_MPI
    , comm
#endif
    );
    dg::x::DVec vol3d = dg::create::volume( g3d);
    //create magnetic field
    const dg::geo::TokamakMagneticField mag = dg::geo::createCircularField( R_0, I_0);
    auto bhat = dg::geo::createBHat(mag);
    auto ff = dg::geo::TestFunctionDirNeu(mag);
    const dg::x::DVec fun = dg::pullback( ff, g3d);
    dg::x::DVec derivative(fun);
    const dg::x::DVec divb = dg::pullback( dg::geo::Divb(mag), g3d);
    const dg::x::DVec sol0 = dg::pullback( dg::geo::DsFunction<dg::geo::TestFunctionDirNeu>(mag,ff), g3d);
    const dg::x::DVec sol1 = dg::pullback( dg::geo::DssFunction<dg::geo::TestFunctionDirNeu>(mag,ff), g3d);
    const dg::x::DVec sol2 = dg::pullback( dg::geo::DsDivFunction<dg::geo::TestFunctionDirNeu>(mag,ff), g3d);
    const dg::x::DVec sol3 = dg::pullback( dg::geo::DsDivDsFunction<dg::geo::TestFunctionDirNeu>(mag,ff), g3d);
    const dg::x::DVec sol4 = dg::pullback( dg::geo::OMDsDivDsFunction<dg::geo::TestFunctionDirNeu>(mag,ff), g3d);
    std::vector<std::tuple<std::string, std::array<const dg::x::DVec*,2>, double>> names{
         {"forward",{&fun,&sol0}, 0.21},          {"backward",{&fun,&sol0}, 0.21},
         {"forward2",{&fun,&sol0}, 0.17},         {"backward2",{&fun,&sol0}, 0.17},
         {"centered",{&fun,&sol0}, 0.08},         {"centered_bc_along",{&fun,&sol0}, 0.08},
         {"dss",{&fun,&sol1}, 0.36},              {"dss_bc_along",{&fun,&sol1}, 0.24},
         {"divForward",{&fun,&sol2}, 0.21},       {"divBackward",{&fun,&sol2}, 0.21},
         {"divCentered",{&fun,&sol2}, 0.08},      {"directLap",{&fun,&sol3}, 0.36},
         {"directLap_bc_along",{&fun,&sol3}, 0.24}//, {"invCenteredLap",{&sol4,&fun}} //inversion may take too long for test
    };
    ///##########################################################///
    //create Fieldaligned object and construct DS from it
    dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>  dsFA(
            bhat, g3d, bcx, bcx, dg::geo::NoLimiter(), 1e-8, mx[0], mx[1],
            -1, method, false);
    dg::geo::DS<dg::x::aProductGeometry3d, dg::x::IDMatrix, dg::x::DVec>
        ds( dsFA );
    for( const auto& tuple :  names)
    {
        std::string name = std::get<0>(tuple);
        const dg::x::DVec& function = *std::get<1>(tuple)[0];
        const dg::x::DVec& solution = *std::get<1>(tuple)[1];
        if( name.find("inv") != std::string::npos ||
                name.find( "div") != std::string::npos)
            callDS( ds, name, function, derivative, max_iter,1e-8);
        else
        {
            // test aliasing
            dg::blas1::copy( function, derivative);
            callDS( ds, name, derivative, derivative, max_iter,1e-8);
        }
        double sol = dg::blas2::dot( vol3d, solution);
        dg::blas1::axpby( 1., solution, -1., derivative);
        double norm = dg::blas2::dot( derivative, vol3d, derivative);
        INFO("    "<<name<<":" <<std::setw(18-name.size())
                  <<" "<<sqrt(norm/sol));
        CHECK( sqrt( norm/sol) <= std::get<2>(tuple));
    }
}

TEST_CASE("Fieldaligned")
{
#ifdef WITH_MPI
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {false,false,true});
#endif
    const dg::x::CylindricalGrid3d g3d( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        3, 26, 26, 20, dg::NEU, dg::NEU, dg::PER
#ifdef WITH_MPI
    , comm
#endif
    );
    const dg::geo::TokamakMagneticField mag = dg::geo::createCircularField( R_0, I_0);
    auto bhat = dg::geo::createBHat(mag);
    std::string method = "dg";
    dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>  dsFA(
            bhat, g3d, dg::NEU, dg::NEU, dg::geo::NoLimiter(), 1e-8, 10, 10,
            -1, method, false);
    dg::geo::DS<dg::x::aProductGeometry3d, dg::x::IDMatrix, dg::x::DVec>
        ds( dsFA );
    dg::x::DVec derivative = dg::evaluate( dg::zero, g3d);
    dg::x::DVec vol3d = dg::create::volume( g3d);
    ///##########################################################///
    INFO( "FIELDALIGNED EVALUATION of a Gaussian");
    dg::Gaussian init0(R_0+0.5, 0, 0.2, 0.2, 1);
    dg::GaussianZ modulate(0., M_PI/3., 1);
    dg::x::DVec aligned = dg::geo::fieldaligned_evaluate( g3d, bhat, init0, modulate, 10, 2);
    ds.ds( dg::centered, aligned, derivative);
    double norm = dg::blas2::dot(vol3d, derivative);
    INFO( "# Norm Centered Derivative "<<sqrt( norm)<<" (compare with that of ds_mpit)");
    CHECK( norm < 0.099);
    aligned = dsFA.evaluate( init0, modulate, 10, 2);
    ds.ds( dg::centered, aligned, derivative);
    norm = dg::blas2::dot(vol3d, derivative);
    INFO( "# Norm Centered Derivative "<<sqrt( norm)<<" (compare with that of ds_mpit)");
    CHECK( norm < 0.099);
    ///##########################################################///
    INFO( "STAGGERED GRID DERIVATIVE");
    auto ff = dg::geo::TestFunctionDirNeu(mag);
    const dg::x::DVec fun = dg::pullback( ff, g3d);
    const dg::x::DVec sol0 = dg::pullback( dg::geo::DsFunction<dg::geo::TestFunctionDirNeu>(mag,ff), g3d);
    dg::x::DVec zMinus(fun), eMinus(fun), zPlus(fun), ePlus(fun), eZero(fun);
    dg::x::DVec funST(fun);
    dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>  dsFAST(
            bhat, g3d, dg::NEU, dg::NEU, dg::geo::NoLimiter(), 1e-8, 10, 10,
            g3d.hz()/2., method, false);

    dsFAST( dg::geo::zeroMinus, fun, zMinus);
    dsFAST( dg::geo::einsPlus,  fun, ePlus);
    dg::geo::assign_bc_along_field_1st( dsFAST, zMinus, ePlus, zMinus, ePlus,
        dg::NEU, {0,0});
    dg::geo::ds_average( dsFAST, 1., zMinus, ePlus, 0., funST);
    dsFAST( dg::geo::zeroPlus, funST, zPlus);
    dsFAST( dg::geo::einsMinus, funST, eMinus);
    dg::geo::assign_bc_along_field_1st( dsFAST, eMinus, zPlus, eMinus, zPlus,
        dg::NEU, {0,0});
    dg::geo::ds_slope( dsFAST, 1., eMinus, zPlus, 0., derivative);
    double sol = dg::blas2::dot( vol3d, sol0);
    dg::blas1::axpby( 1., sol0, -1., derivative);
    norm = dg::blas2::dot( derivative, vol3d, derivative);
    std::string name = "forward";
    INFO("    "<<name<<":" <<std::setw(18-name.size())
              <<" "<<sqrt(norm/sol));
    CHECK( sqrt(norm/sol) < 0.0711);

    // now try the adjoint direction (should be exactly the same result)
    dsFAST( dg::geo::zeroPlus, fun, zPlus);
    dsFAST( dg::geo::einsMinus, fun, eMinus);
    dg::geo::assign_bc_along_field_1st( dsFAST, eMinus, zPlus, eMinus, zPlus,
        dg::NEU, {0,0});
    dg::geo::ds_average( dsFAST, 1., eMinus, zPlus, 0., funST);
    dsFAST( dg::geo::einsPlus, funST, ePlus);
    dsFAST( dg::geo::zeroMinus, funST, zMinus);
    dg::geo::assign_bc_along_field_1st( dsFAST, zMinus, ePlus, zMinus, ePlus,
        dg::NEU, {0,0});
    dg::geo::ds_slope( dsFAST, 1., zMinus, ePlus, 0., derivative);
    dg::blas1::axpby( 1., sol0, -1., derivative);
    norm = dg::blas2::dot( derivative, vol3d, derivative);
    name = "backward";
    INFO("    "<<name<<":" <<std::setw(18-name.size())
              <<" "<<sqrt(norm/sol));
    CHECK( sqrt(norm/sol) < 0.0711);
}

// The member functions are documented to do the same as the freestanding ones
TEST_CASE( "DS boundary conditions along fieldline")
{
#ifdef WITH_MPI
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {false,false,true});
#endif
    const dg::x::CylindricalGrid3d g3d( R_0-a, R_0+a, -a, a, 0, 2.*M_PI,
        3, 8, 8, 4, dg::NEU, dg::NEU, dg::PER
#ifdef WITH_MPI
    , comm
#endif
    );
    const dg::geo::TokamakMagneticField mag = dg::geo::createCircularField( R_0, I_0);
    const std::string method = GENERATE( "dg", "linear-nearest");
    const dg::bc bound = GENERATE( dg::DIR, dg::NEU);
    INFO( "method "<<method<<" bound "<<dg::bc2str(bound));
    using FA = dg::geo::Fieldaligned<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec>;
    FA fa( dg::geo::createBHat(mag), g3d, dg::NEU, dg::NEU, dg::geo::NoLimiter(), 1e-8, 3, 3, -1, method, false);
    dg::geo::DS<dg::x::aProductGeometry3d,dg::x::IDMatrix,dg::x::DVec> ds( fa);
    const dg::x::DVec fun = dg::evaluate( []( double R, double Z, double phi){
        return exp(R-R_0)*(1+Z*Z)*(1+sin(phi));}, g3d);
    const std::array<double,2> value{1.,2.};
    dg::x::DVec fm( fun), fo( fun), fp( fun), member( fun), free( fun);
    fa( dg::geo::einsMinus, fun, fm);
    fa( dg::geo::zeroForw, fun, fo);
    fa( dg::geo::einsPlus, fun, fp);
    dg::geo::assign_bc_along_field_2nd( fa, fm, fo, fp, fm, fp, bound, value);
    auto difference = [&](){
        dg::blas1::axpby( 1., member, -1., free);
        return dg::blas1::dot( free, free);
    };
    SECTION( "centered")
    {
        ds.centered_bc_along_field( 1., fun, 0., member, bound, value);
        dg::geo::ds_centered( fa, 1., fm, fp, 0., free);
        CHECK( difference() == 0);
    }
    SECTION( "dss")
    {
        ds.dss_bc_along_field( 1., fun, 0., member, bound, value);
        dg::geo::dss_centered( fa, 1., fm, fo, fp, 0., free);
        CHECK( difference() == 0);
    }
    SECTION( "dssd")
    {
        ds.dssd_bc_along_field( 1., fun, 0., member, bound, value);
        dg::geo::dssd_centered( fa, 1., fm, fo, fp, 0., free);
        CHECK( difference() == 0);
    }
}

// Provides what the assign_bc functions need from Fieldaligned for a single point
struct SinglePoint
{
    SinglePoint( double delta, double hm, double hp, double bm, double bo, double bp) :
        m_delta(delta), m_hbm(1,hm), m_hbp(1,hp), m_bbm(1,bm), m_bbo(1,bo), m_bbp(1,bp){}
    double deltaPhi() const{ return m_delta;}
    const dg::HVec& hbm() const{ return m_hbm;}
    const dg::HVec& hbp() const{ return m_hbp;}
    const dg::HVec& bbm() const{ return m_bbm;}
    const dg::HVec& bbo() const{ return m_bbo;}
    const dg::HVec& bbp() const{ return m_bbp;}
    private:
    double m_delta;
    dg::HVec m_hbm, m_hbp, m_bbm, m_bbo, m_bbp;
};

TEST_CASE( "Boundary conditions along fieldline")
{
    // s=0 is the center, the boundaries lie at s=-hbm and s=hbp and the minus
    // and plus points at s=-delta and s=+delta (they are outside if the mask says so)
    const double delta = 0.3, hbm = 0.2, hbp = 0.1, eps = 1e-14;
    const double c0 = 1.3, c1 = 0.7, c2 = -1.1;
    auto poly = [=]( double s){ return c0 + c1*s + c2*s*s;};
    auto dpoly = [=]( double s){ return c1 + 2.*c2*s;};
    auto line = [=]( double s){ return c0 + c1*s;};
    enum Mask{ minus, both, plus};
    const Mask mask = GENERATE( minus, both, plus);
    INFO( "Mask "<<(mask == minus ? "minus" : mask == both ? "both" : "plus"));
    const SinglePoint fa( delta, hbm, hbp, mask == minus, mask == both, mask == plus);
    dg::HVec fmg(1), fpg(1);

    SECTION( "2nd order reproduces a parabola")
    {
        const dg::HVec fm(1,poly(-delta)), fo(1,poly(0)), fp(1,poly(delta));
        dg::geo::assign_bc_along_field_2nd( fa, fm, fo, fp, fmg, fpg, dg::DIR,
                {poly(-hbm), poly(hbp)});
        CHECK_THAT( fmg[0], Catch::Matchers::WithinAbs( poly(-delta), eps));
        CHECK_THAT( fpg[0], Catch::Matchers::WithinAbs( poly(+delta), eps));
        dg::geo::assign_bc_along_field_2nd( fa, fm, fo, fp, fmg, fpg, dg::NEU,
                {dpoly(-hbm), dpoly(hbp)});
        CHECK_THAT( fmg[0], Catch::Matchers::WithinAbs( poly(-delta), eps));
        CHECK_THAT( fpg[0], Catch::Matchers::WithinAbs( poly(+delta), eps));
    }
    SECTION( "1st order reproduces a line")
    {
        const dg::HVec fm(1,line(-delta)), fp(1,line(delta));
        dg::geo::assign_bc_along_field_1st( fa, fm, fp, fmg, fpg, dg::DIR,
                {line(-hbm), line(hbp)});
        CHECK_THAT( fmg[0], Catch::Matchers::WithinAbs( line(-delta), eps));
        CHECK_THAT( fpg[0], Catch::Matchers::WithinAbs( line(+delta), eps));
        dg::geo::assign_bc_along_field_1st( fa, fm, fp, fmg, fpg, dg::NEU, {c1, c1});
        CHECK_THAT( fmg[0], Catch::Matchers::WithinAbs( line(-delta), eps));
        CHECK_THAT( fpg[0], Catch::Matchers::WithinAbs( line(+delta), eps));
    }
    SECTION( "Only points outside are replaced")
    {
        // These input values are inconsistent with the boundary values
        const dg::HVec fm(1,5.), fp(1,7.);
        dg::geo::assign_bc_along_field_1st( fa, fm, fp, fmg, fpg, dg::DIR, {1., 2.});
        CHECK( (fmg[0] == 5.) == (mask == plus));
        CHECK( (fpg[0] == 7.) == (mask == minus));
        dg::geo::assign_bc_along_field_1st( fa, fm, fp, fmg, fpg, dg::NEU, {1., 2.});
        CHECK( (fmg[0] == 5.) == (mask != minus));
        CHECK( (fpg[0] == 7.) == (mask != plus));
    }
    SECTION( "Swap perpendicular boundary condition")
    {
        const dg::HVec fm(1,5.), fp(1,7.);
        dg::geo::swap_bc_perp( fa, fm, fp, fmg, fpg);
        CHECK( fmg[0] == (mask == plus ? +5. : -5.));
        CHECK( fpg[0] == (mask == minus ? +7. : -7.));
    }
}
