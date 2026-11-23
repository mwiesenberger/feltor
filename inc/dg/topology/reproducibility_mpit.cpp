#include <iostream>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <limits>
#include <mpi.h>
#include "dg/backend/mpi_init.h"
#include "dg/blas.h"
#include "mpi_projection.h"
#include "mpi_evaluation.h"
#include "average.h"

#include "catch2/catch_all.hpp"

// The MPI matrices must give bitwise the same result as the shared memory
// matrices of the global grid on the same backend, for any number of ranks
// and any domain decomposition.

static double test_function( double x, double y)
{
    // the small high frequency term makes the summation order visible in the
    // last bits
    return sin(M_PI*x/2.)*sin(M_PI*y/2.) + 1e-3*cos( 13.*x*y) + 1.;
}

// Number of representable doubles between a and b (0 if bitwise equal)
static std::uint64_t ulp_distance( double a, double b)
{
    std::int64_t ia, ib;
    std::memcpy( &ia, &a, sizeof(double));
    std::memcpy( &ib, &b, sizeof(double));
    // map the sign magnitude representation to a monotonic integer scale
    if( ia < 0)
        ia = std::numeric_limits<std::int64_t>::min() - ia;
    if( ib < 0)
        ib = std::numeric_limits<std::int64_t>::min() - ib;
    return ia > ib ? std::uint64_t(ia) - std::uint64_t(ib)
                   : std::uint64_t(ib) - std::uint64_t(ia);
}

// Compare the local MPI result with the expected local values: they may
// differ by at most max_ulp ulp (0 means bitwise equal)
static void check_ulp( const dg::HVec& local, const dg::HVec& expected,
    unsigned max_ulp)
{
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    unsigned mismatches = 0;
    std::uint64_t worst_ulp = 0;
    unsigned worst = 0;
    for( unsigned i=0; i<local.size(); i++)
    {
        std::uint64_t dist = ulp_distance( local[i], expected[i]);
        if( dist > 0)
            mismatches++;
        if( dist > worst_ulp)
        {
            worst_ulp = dist;
            worst = i;
        }
    }
    INFO( "Rank "<<rank<<": "<<mismatches<<" of "<<local.size()
        <<" values differ, at most by "<<worst_ulp<<" ulp (allowed "<<max_ulp
        <<"), at local index "<<worst);
    if( local.size() > 0)
        CHECK_THAT( local[worst], Catch::Matchers::WithinULP( expected[worst],
            max_ulp));
}


// y = M x on the MPI grids and on the global grids, with the MPI matrix and
// the global matrix copied to the device, and compare bitwise.
static void check_grid_matrix( const dg::MIHMatrix& mpi_host,
    const dg::IHMatrix& global_host, const dg::MPIGrid2d& g_new,
    const dg::MPIGrid2d& g_old)
{
    dg::MIDMatrix mpi_mat = mpi_host;
    dg::IDMatrix global_mat = global_host;

    dg::HVec global_x = dg::evaluate( test_function, g_old.global());
    dg::MDVec x = dg::global2local( global_x, g_old);
    dg::MDVec y = dg::construct<dg::MDVec>( dg::evaluate( dg::zero, g_new));
    dg::blas2::symv( mpi_mat, x, y);

    dg::DVec gx = global_x;
    dg::DVec gy( g_new.global().size());
    dg::blas2::symv( global_mat, gx, gy);

    check_ulp( y.data(), dg::global2local( dg::HVec( gy), g_new).data(), 0);
}

TEST_CASE( "MPI matrices agree with shared memory matrices")
{
    int size;
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {1,0});
    INFO( "Number of ranks "<<size);

    // The cell numbers are not multiples of the number of ranks, so that
    // the local grids have different sizes
    dg::MPIGrid2d g_coarse( 0.1, 1.7, -0.3, 1.2, 3, 7, 11, dg::PER, dg::DIR, comm);

    SECTION( "Interpolation from a coarse to a finer grid")
    {
        std::string method = GENERATE( "dg", "linear", "cubic", "nearest");
        INFO( "Method "<<method);
        dg::MPIGrid2d g_fine( 0.1, 1.7, -0.3, 1.2, 4, 13, 9, dg::PER, dg::DIR, comm);
        auto mpi_host = dg::create::interpolation( g_fine, g_coarse, method);
        auto global_host = dg::create::interpolation( g_fine.global(),
                g_coarse.global(), method);
        check_grid_matrix( mpi_host, global_host, g_fine, g_coarse);
    }
    SECTION( "Projection from a fine to a coarser grid")
    {
        dg::MPIGrid2d g_fine = g_coarse;
        g_fine.multiplyCellNumbers( 3, 2);
        auto mpi_host = dg::create::projection( g_coarse, g_fine);
        auto global_host = dg::create::projection( g_coarse.global(),
                g_fine.global());
        check_grid_matrix( mpi_host, global_host, g_coarse, g_fine);
    }
    SECTION( "Interpolation to given points")
    {
        // Every rank has its own points; the shift makes the interpolation
        // communicate
        std::string method = GENERATE( "dg", "linear", "cubic", "nearest");
        INFO( "Method "<<method);
        dg::HVec xs = dg::evaluate( dg::cooX2d, g_coarse).data();
        dg::HVec ys = dg::evaluate( dg::cooY2d, g_coarse).data();
        for( unsigned i=0; i<xs.size(); i++)
        {
            xs[i] = 0.1 + fmod( xs[i] - 0.1 + 0.37, 1.6);
            ys[i] = -0.3 + 0.97*(ys[i]+0.3);
        }
        auto mpi_host = dg::create::interpolation( xs, ys, g_coarse,
                dg::PER, dg::DIR, method);
        auto global_host = dg::create::interpolation( xs, ys,
                g_coarse.global(), dg::PER, dg::DIR, method);
        dg::MIDMatrix mpi_mat = mpi_host;
        dg::IDMatrix global_mat = global_host;

        dg::HVec global_x = dg::evaluate( test_function, g_coarse.global());
        dg::MDVec x = dg::global2local( global_x, g_coarse);
        dg::MDVec y( dg::DVec( xs.size()), g_coarse.communicator());
        dg::blas2::symv( mpi_mat, x, y);

        dg::DVec gx = global_x;
        dg::DVec gy( xs.size());
        dg::blas2::symv( global_mat, gx, gy);
        // the rows are the same local points
        check_ulp( y.data(), gy, 0);
    }
}

// dg::Average, with and without extension, and dg::create::reduction along
// the given axes on the MPI grid and on the global grid, compared bitwise
template<size_t Nd, size_t Md>
static void check_average( const dg::RealMPIGrid<double,Nd>& g,
    std::array<unsigned,Md> axes, const dg::HVec& global_x)
{
    std::array<dg::MPIGrid1d, Nd-Md> gs;
    unsigned k = 0;
    for( unsigned u=0; u<Nd; u++)
        if( std::find( axes.begin(), axes.end(), u) == axes.end())
            gs[k++] = g.grid(u);
    dg::RealMPIGrid<double,Nd-Md> g_red( gs);
    dg::MDVec x = dg::global2local( global_x, g);
    dg::DVec gx = global_x;

    dg::Average< dg::MIDMatrix, dg::MDVec> mpi_avg( g, axes);
    dg::Average< dg::IDMatrix, dg::DVec> avg( g.global(), axes);
    for( bool extend : { true, false})
    {
        INFO( "Extend "<<extend);
        dg::MDVec y( x);
        dg::DVec gy( gx);
        mpi_avg( x, y, extend);
        avg( gx, gy, extend);
        if( extend)
            check_ulp( y.data(), dg::global2local( dg::HVec( gy), g).data(), 0);
        else
            check_ulp( y.data(), dg::global2local( dg::HVec( gy), g_red).data(), 0);
    }

    dg::MIDMatrix mpi_red = dg::create::reduction( axes, g);
    dg::IDMatrix red = dg::create::reduction( axes, g.global());
    dg::MDVec y = dg::construct<dg::MDVec>( dg::evaluate( dg::zero, g_red));
    dg::DVec gy( g_red.global().size());
    dg::blas2::symv( mpi_red, x, y);
    dg::blas2::symv( red, gx, gy);
    INFO( "Reduction");
    check_ulp( y.data(), dg::global2local( dg::HVec( gy), g_red).data(), 0);
}

TEST_CASE( "MPI Average agrees with shared memory Average")
{
    int size;
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    // At least as many cells as ranks on every axis, so that every remaining
    // axis can hold all processes, and not a multiple of the number of ranks
    const unsigned Nx = 2*size+1, Ny = 2*size+3, Nz = 2*size+5;
    SECTION( "Average in 2d")
    {
        // Average along x or y, with the ranks distributed only along x, only
        // along y, or along both
        auto axis = GENERATE( 0u, 1u);
        std::array<int,2> dims = GENERATE( std::array<int,2>{0,1},
            std::array<int,2>{1,0}, std::array<int,2>{0,0});
        INFO( "Number of ranks "<<size<<", average along axis "<<axis
            <<", dims "<<dims[0]<<" "<<dims[1]);
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD,
            {dims[0], dims[1]}, {1,0});
        dg::MPIGrid2d g( 0.1, 1.7, -0.3, 1.2, 3, Nx, Ny, dg::PER, dg::DIR, comm);
        check_average( g, std::array{axis},
            dg::evaluate( test_function, g.global()));
    }
    SECTION( "Average in 3d")
    {
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {1,0,0});
        dg::MPIGrid3d g( 0.1, 1.7, -0.3, 1.2, 0.2, 0.9, 2, Nx, Ny, Nz, dg::PER,
            dg::DIR, dg::NEU, comm);
        dg::HVec global_x = dg::evaluate( []( double x, double y, double z){
                return test_function( x, y)*(1.+0.1*z);}, g.global());
        auto axis = GENERATE( 0u, 1u, 2u);
        INFO( "Number of ranks "<<size<<", average along axis "<<axis
            <<" and along all axes except "<<axis);
        check_average( g, std::array{axis}, global_x);
        std::array<unsigned,2> pair = axis == 0 ? std::array{1u,2u} :
            ( axis == 1 ? std::array{0u,2u} : std::array{0u,1u});
        check_average( g, pair, global_x);
    }
    SECTION( "More ranks than cells")
    {
        // With more than 3 ranks some processes get none of the 3 cells along
        // x: in the redistribution for the average along y, and with dims
        // {size,1} also in the grid itself
        std::array<int,2> dims = GENERATE( std::array<int,2>{1,0},
            std::array<int,2>{0,1});
        auto axis = GENERATE( 0u, 1u);
        INFO( "Number of ranks "<<size<<", average along axis "<<axis
            <<", dims "<<dims[0]<<" "<<dims[1]);
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD,
            {dims[0], dims[1]}, {1,0});
        dg::MPIGrid2d g( 0.1, 1.7, -0.3, 1.2, 3, 3, 11, dg::PER, dg::DIR, comm);
        check_average( g, std::array{axis},
            dg::evaluate( test_function, g.global()));
    }
}
