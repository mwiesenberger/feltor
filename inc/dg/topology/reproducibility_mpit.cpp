#include <iostream>
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
// and any domain decomposition. The known exception is checked against a
// bound in units in the last place (ulp): dg::Average adds the partial sums
// of the ranks along the averaged directions with MPI_Allreduce

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

TEST_CASE( "MPI Average agrees with shared memory Average")
{
    int size;
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    // Average along x or y, with the ranks distributed only along x, only
    // along y, or along both. If the averaged direction is not distributed
    // the MPI Average needs no reduction across ranks
    auto axis = GENERATE( dg::coo2d::x, dg::coo2d::y);
    std::array<int,2> dims = GENERATE( std::array<int,2>{0,1},
        std::array<int,2>{1,0}, std::array<int,2>{0,0});
    INFO( "Number of ranks "<<size<<", average along "
        <<(axis == dg::coo2d::x ? "x" : "y")<<", dims "<<dims[0]<<" "<<dims[1]);
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD,
        {dims[0], dims[1]}, {1,0});
    dg::MPIGrid2d g( 0.1, 1.7, -0.3, 1.2, 3, 7, 11, dg::PER, dg::DIR, comm);

    dg::HVec global_x = dg::evaluate( test_function, g.global());
    dg::MDVec x = dg::global2local( global_x, g);
    dg::MDVec y( x);
    dg::Average< dg::MIDMatrix, dg::MDVec> mpi_avg( g, axis);
    mpi_avg( x, y, true);

    dg::DVec gx = global_x, gy( gx);
    dg::Average< dg::IDMatrix, dg::DVec> avg( g.global(), axis);
    avg( gx, gy, true);

    // Number of ranks along the averaged direction: if larger than 1 the
    // partial sums of the ranks are added with MPI_Allreduce in an order
    // defined by the MPI library
    int P = 1;
    MPI_Comm_size( g.comm( axis == dg::coo2d::x ? 0 : 1), &P);
    unsigned max_ulp = P == 1 ? 0 : 4*P;
    INFO( "Ranks along the averaged direction "<<P);
    check_ulp( y.data(), dg::global2local( dg::HVec( gy), g).data(), max_ulp);
}
