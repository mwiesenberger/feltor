#include <iostream>
#ifdef WITH_MPI
#include <mpi.h>
#include "../backend/mpi_init.h"
#include "mpi_evaluation.h"
#endif
#include "evaluation.h"
#include "split_and_join.h"
#include "dg/blas.h"

#include "catch2/catch_all.hpp"

static double function( double x, double y, double z) { return x + 10.*y + 100.*z;}

TEST_CASE( "Split a 3d vector into planes")
{
#ifdef WITH_MPI
    MPI_Comm comm3d = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {0,0,0});
#endif
    // Nx, Ny, Nz divisible by the process numbers for 1 to 8 ranks
    const dg::x::Grid3d g( 0, 1, 0, 2, 0, 3, 2, 6, 6, 8
#ifdef WITH_MPI
        , comm3d
#endif
    );
#ifdef WITH_MPI
    using View = dg::get_mpi_view_type<dg::x::DVec>;
    const dg::x::Grid2d g2d( 0, 1, 0, 2, 2, 6, 6, g.get_perp_comm());
    const dg::Grid3d local = g.local();
#else
    using View = dg::View<dg::x::DVec>;
    const dg::x::Grid2d g2d( 0, 1, 0, 2, 2, 6, 6);
    const dg::Grid3d local = g;
#endif
    const thrust::host_vector<double> zs = local.abscissas(2);
    dg::x::DVec vec = dg::evaluate( function, g);
    const dg::x::DVec copy = vec;

    // Plane k holds the function on the k-th (local) z-plane and writing to
    // the planes writes to vec
    auto check_planes = [&]( std::vector<View>& planes)
    {
        REQUIRE( planes.size() == local.shape(2));
        for( unsigned k=0; k<planes.size(); k++)
        {
            double z = zs[k];
            dg::x::DVec diff = dg::evaluate( [z]( double x, double y){
                    return function( x, y, z);}, g2d);
            dg::blas1::axpby( 1., planes[k], -1., diff);
            // The compiler may round differently in 2d and 3d (fma contraction)
            INFO( "Plane "<<k);
            CHECK( sqrt( dg::blas1::dot( diff, diff)) < 1e-12);
            dg::blas1::scal( planes[k], 2.);
        }
        dg::blas1::axpby( 2., copy, -1., vec);
        CHECK( dg::blas1::dot( vec, vec) == 0);
    };
    SECTION( "Construct version")
    {
        std::vector<View> planes = dg::split( vec, g);
#ifdef WITH_MPI
        for( auto& plane : planes)
        {
            int result;
            MPI_Comm_compare( plane.communicator(), g.get_perp_comm(), &result);
            CHECK( result == MPI_IDENT);
        }
#endif
        check_planes( planes);
    }
    SECTION( "Fast version")
    {
        std::vector<View> planes( local.shape(2));
#ifdef WITH_MPI
        for( auto& plane : planes)
            plane.set_communicator( g.get_perp_comm());
#endif
        dg::split( vec, planes, g);
        check_planes( planes);
    }
}
