#include <iostream>
#ifdef WITH_MPI
#include <mpi.h>
#include "../backend/mpi_init.h"
#include "mpi_evaluation.h"
#include "mpi_derivatives.h"
#include "mpi_weights.h"
#endif

#include "dg/blas.h"
#include "dx.h"
#include "derivatives.h"
#include "derivativesT.h"
#include "evaluation.h"
#include "weights.h"

#include "catch2/catch_all.hpp"

static double function( double x) { return sin(x);}
static double derivative( double x) { return cos(x);}

using Vector = dg::x::HVec;
using Matrix = dg::x::HMatrix;

TEST_CASE( "Dx")
{
    unsigned n=7, N=33;
#ifdef WITH_MPI
    MPI_Comm comm1d = dg::mpi_cart_create( MPI_COMM_WORLD, {0}, {0});
    MPI_Comm comm1dPER = dg::mpi_cart_create( MPI_COMM_WORLD, {0}, {1});
#endif
    auto bcx = GENERATE( dg::PER, dg::DIR, dg::NEU, dg::DIR_NEU, dg::NEU_DIR);
    auto dir = GENERATE( dg::centered, dg::forward, dg::backward);
    dg::x::Grid1d g1d( 0.1, 2*M_PI+0.1, n, N, bcx
#ifdef WITH_MPI
            , bcx == dg::PER ? comm1dPER : comm1d
#endif
    );

    if( bcx == dg::DIR)
        g1d.set_pq( {0}, {M_PI});
    else if( bcx == dg::NEU)
        g1d.set_pq( {M_PI/2.}, {3.*M_PI/2.});
    else if( bcx == dg::DIR_NEU)
        g1d.set_pq( {0.}, {M_PI/2.});
    else if( bcx == dg::NEU_DIR)
        g1d.set_pq( {M_PI/2.}, {M_PI});

    const Vector func = dg::evaluate( function, g1d);
    const Vector w1d = dg::create::weights( g1d);
    const Vector deri = dg::evaluate( derivative, g1d);
    const Vector null = dg::evaluate( dg::zero, g1d);
    Matrix js = dg::create::jumpX( g1d, g1d.bcx());
    Vector error = func;
    dg::blas2::symv( js, func, error);
    dg::blas1::axpby( 1., null , -1., error);
    double dist = sqrt( dg::blas2::dot( w1d, error));
    INFO("Distance to true solution (jump     ): "<<dist);
    CHECK( dist < 1e-10);

    Matrix hs = dg::create::dx( g1d, g1d.bcx(), dir);

    dg::blas2::symv( hs, func, error);
    dg::blas1::axpby( 1., deri, -1., error);
    dist = sqrt( dg::blas2::dot( w1d, error));
    INFO( "Distance to true solution "<<dg::direction2str(dir )<<dist);
    CHECK( dist < 1e-10);
    //for periodic bc | dirichlet bc
    //n = 1 -> p = 2      2
    //n = 2 -> p = 1      1
    //n = 3 -> p = 3      3
    //n = 4 -> p = 3      3
    //n = 5 -> p = 5      5
}

static double quadratic( double x, dg::bc bcx)
{
    // satisfies the Dirichlet conditions on [0,1] where they apply
    switch( bcx)
    {
        case dg::NEU: return (x+1.)*(x+1.);
        case dg::DIR_NEU: return x*(x+1.);
        case dg::NEU_DIR: return x*x-1.;
        default: return x*(1.-x);
    }
}
static double dquadratic( double x, dg::bc bcx)
{
    switch( bcx)
    {
        case dg::NEU: return 2.*(x+1.);
        case dg::DIR_NEU: return 2.*x+1.;
        case dg::NEU_DIR: return 2.*x;
        default: return 1.-2.*x;
    }
}

TEST_CASE( "Dx with a single cell")
{
    // With n=3 quadratics are represented exactly and the discrete
    // derivative of a function compatible with the boundary conditions is
    // exact also with a single cell
    auto bcx = GENERATE( dg::PER, dg::DIR, dg::NEU, dg::DIR_NEU, dg::NEU_DIR);
    auto dir = GENERATE( dg::centered, dg::forward, dg::backward);
    auto N = GENERATE( 1u, 3u);
    dg::Grid1d g1d( 0., 1., 3, N, bcx);
    const dg::HVec func = dg::evaluate( [bcx](double x){ return
            quadratic(x, bcx);}, g1d);
    const dg::HVec deri = dg::evaluate( [bcx](double x){ return
            dquadratic(x, bcx);}, g1d);
    dg::HVec result( func);
    dg::HMatrix jump = dg::create::jumpX( g1d, bcx);
    dg::blas2::symv( jump, func, result);
    INFO( "N "<<N<<" bc "<<dg::bc2str(bcx)<<" "<<dg::direction2str(dir));
    for( unsigned i=0; i<result.size(); i++)
        CHECK( fabs( result[i]) < 1e-13);
    dg::HMatrix dx = dg::create::dx( g1d, bcx, dir);
    dg::blas2::symv( dx, func, result);
    for( unsigned i=0; i<result.size(); i++)
        CHECK( fabs( result[i] - deri[i]) < 1e-13);
}
