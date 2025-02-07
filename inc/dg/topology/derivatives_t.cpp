#include <iostream>
#include <cmath>
#ifdef WITH_MPI
#include <mpi.h>
#include "../backend/mpi_init.h"
#include "mpi_evaluation.h"
#include "mpi_derivatives.h"
#include "mpi_weights.h"
#endif
#include "dg/blas.h"
#include "dg/functors.h"
#include "evaluation.h"
#include "derivatives.h"
#include "derivativesT.h"

#include "catch2/catch_all.hpp"

using Matrix = dg::x::DMatrix;
using Vector = dg::x::DVec;
using value_t = double;

static value_t zero( value_t x, value_t y) { return 0;}
static value_t sine( value_t x, value_t y) { return sin(x)*sin(y);}
static value_t cosx( value_t x, value_t y) { return cos(x)*sin(y);}
static value_t cosy( value_t x, value_t y) { return cos(y)*sin(x);}
static value_t zero( value_t x, value_t y, value_t z) { return 0;}
static value_t sine( value_t x, value_t y, value_t z) { return sin(x)*sin(y)*sin(z);}
static value_t cosx( value_t x, value_t y, value_t z) { return cos(x)*sin(y)*sin(z);}
static value_t cosy( value_t x, value_t y, value_t z) { return cos(y)*sin(x)*sin(z);}
static value_t cosz( value_t x, value_t y, value_t z) { return cos(z)*sin(x)*sin(y);}

// It seems going from g++-11 to g++-13 changes the implementation of sin and cos

TEST_CASE( "Derivatives")
{
#ifdef WITH_MPI
    int size;
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    std::vector<int> dims = {0,0,0};
    MPI_Dims_create( size, 3, &dims[0]);
    auto i = GENERATE( 0,1,2,3,4,5);
    std::sort( dims.begin(), dims.end());
    for( int u=0; u<i; u++)
        std::next_permutation( dims.begin(), dims.end());
    INFO( "Permutation of dims "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]);
    std::vector<int> dims2d = {dims[0], dims[1]};
    MPI_Comm comm3d = dg::mpi_cart_create( MPI_COMM_WORLD, dims, {0, 1, 0});
    auto comms1d = dg::mpi_cart_split( comm3d);
    MPI_Comm comm2d = dg::mpi_cart_kron( {comms1d[0], comms1d[1]});
#endif
    INFO("This program tests the creation and application of two-dimensional "
            <<"and three-dimensional derivatives!");
    unsigned n = 3, Nx = 24, Ny = 28, Nz = 100;
    dg::bc bcx=dg::DIR, bcy=dg::PER, bcz=dg::NEU_DIR;
    SECTION( "Two dimensional")
    {
        INFO("On Grid "<<n<<" x "<<Nx<<" x "<<Ny);
        dg::x::RealGrid2d<value_t> g2d( 0, M_PI, 0.1, 2*M_PI+0.1, n, Nx, Ny,
                bcx, bcy
#ifdef WITH_MPI
                , comm2d
#endif
                );
        const Vector w2d = dg::create::weights( g2d);

        Matrix dx2 = dg::create::dx( g2d, g2d.bcx(), dg::forward);
        Matrix dy2 = dg::create::dy( g2d, g2d.bcy(), dg::centered);
        Matrix jx2 = dg::create::jumpX( g2d, g2d.bcx());
        Matrix jy2 = dg::create::jumpY( g2d, g2d.bcy());
        Matrix m2[] = {dx2, dy2, jx2, jy2};
        const Vector f2d = dg::evaluate( sine, g2d);
        const Vector dx2d = dg::evaluate( cosx, g2d);
        const Vector dy2d = dg::evaluate( cosy, g2d);
        const Vector null2 = dg::evaluate( zero, g2d);
        Vector sol2[] = {dx2d, dy2d, null2, null2};
#if __GNUC__ >= 13
        int64_t binary2[] = {4562611930300282861,4553674328256673277,4567083257206217158,4574111364446550181};
#else
        int64_t binary2[] = {4562611930300281864,4553674328256556132,4567083257206218817,4574111364446550002};
#endif

        dg::exblas::udouble res;
        INFO("TEST 2D: DX, DY, JX, JY");
        auto i = GENERATE( 0,1,2,3);
        Vector error = sol2[i];
        dg::blas2::symv( -1., m2[i], f2d, 1., error);
        dg::blas1::pointwiseDot( error, error, error);
        value_t norm = sqrt(dg::blas1::dot( w2d, error)); res.d = norm;
        INFO("Distance "<<i<<" to true solution: "<<norm<<"\t"<<res.i<<"\t"<<binary2[i]);
        CHECK( std::abs(res.i - binary2[i]) < 2);
    }
    SECTION( "Three dimensional")
    {
        INFO("On Grid "<<n<<" x "<<Nx<<" x "<<Ny<<" x "<<Nz);
        dg::x::RealGrid3d<value_t> g3d( 0,M_PI, 0.1, 2.*M_PI+0.1, M_PI/2.,M_PI,
                n, Nx, Ny, Nz, bcx, bcy, bcz
#ifdef WITH_MPI
                , comm3d
#endif
                );
        const Vector w3d = dg::create::weights( g3d);
        Matrix dx3 = dg::create::dx( g3d, g3d.bcx(), dg::forward);
        Matrix dy3 = dg::create::dy( g3d, g3d.bcy(), dg::centered);
        Matrix dz3 = dg::create::dz( g3d, g3d.bcz(), dg::backward);
        Matrix jx3 = dg::create::jumpX( g3d, g3d.bcx());
        Matrix jy3 = dg::create::jumpY( g3d, g3d.bcy());
        Matrix jz3 = dg::create::jumpZ( g3d, g3d.bcz());
        Matrix m3[] = {dx3, dy3, dz3, jx3, jy3, jz3};
        const Vector f3d = dg::evaluate( sine, g3d);
        const Vector dx3d = dg::evaluate( cosx, g3d);
        const Vector dy3d = dg::evaluate( cosy, g3d);
        const Vector dz3d = dg::evaluate( cosz, g3d);
        const Vector null3 = dg::evaluate( zero, g3d);
        Vector sol3[] = {dx3d, dy3d, dz3d, null3, null3, null3};
#if __GNUC__ >= 13
        int64_t binary3[] = {4561946736820640320,4553062895410783769,4594213495911299616,4566393134538622288,4573262464593641524,4594304523193682043};
#else
        int64_t binary3[] = {4561946736820639666,4553062895410573431,4594213495911299616,4566393134538626348,4573262464593641240,4594304523193682043};
#endif

        INFO("TEST 3D: DX, DY, DZ, JX, JY, JZ");
        auto i = GENERATE( 0,1,2,3,4,5);
        dg::exblas::udouble res;
        Vector error = sol3[i];
        dg::blas2::symv( -1., m3[i], f3d, 1., error);
        value_t norm = sqrt(dg::blas2::dot( error, w3d, error)); res.d = norm;
        INFO("Distance "<<i<<" to true solution: "<<norm<<"\t"<<res.i<<"\t"<<binary3[i]);
        CHECK( std::abs(res.i - binary3[i]) < 2);
    }
    SECTION( "Symv captures NaN")
    {
        dg::x::RealGrid3d<value_t> g3d( 0,M_PI, 0.1, 2.*M_PI+0.1, M_PI/2.,M_PI,
                n, Nx, Ny, Nz, bcx, bcy, bcz
#ifdef WITH_MPI
                , comm3d
#endif
                );
        Matrix dx3 = dg::create::dx( g3d, g3d.bcx(), dg::forward);
        const Vector f3d = dg::evaluate( sine, g3d);
        Vector error = dg::evaluate( zero, g3d);
#ifdef WITH_MPI
        error.data()[0]  = NAN;
        error.data()[1] = NAN;
#else
        error[0]  = NAN;
        error[10] = NAN;
#endif
        dg::blas2::symv(  dx3, f3d, error);
        dg::x::HVec x( error);
        bool hasnan = dg::blas1::reduce( x, false,
                thrust::logical_or<bool>(), dg::ISNFINITE<double>());
        INFO("Symv contains NaN: "<<std::boolalpha<<hasnan<<" (false)");
        CHECK( not hasnan);
    }
}
