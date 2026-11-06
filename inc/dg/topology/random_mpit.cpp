#include <iostream>
#include <mpi.h>

#include "dg/backend/mpi_init.h"
#include "mpi_evaluation.h"
#include "random.h"

#include "catch2/catch_all.hpp"

TEST_CASE( "MPI uniform real distribution")
{
    // The MPI vector is the distributed version of the shared memory vector
    // for any number of processes
    double a = -1, b = 2;
    SECTION( "2d")
    {
        MPI_Comm comm2d = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {0,1});
        dg::MPIGrid2d g2d( 0, 1, 0, 1, 3, 12, 20, comm2d);
        auto vec = dg::create::uniform_real_distribution( a, b, g2d);
        auto global = dg::create::uniform_real_distribution( a, b, g2d.global());
        auto local = dg::global2local( global, g2d);
        REQUIRE( vec.data().size() == local.data().size());
        for( unsigned u=0; u<vec.data().size(); u++)
        {
            INFO( "Local index "<<u);
            CHECK( vec.data()[u] == local.data()[u]);
        }
    }
    SECTION( "3d")
    {
        MPI_Comm comm3d = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {0,0,1});
        dg::MPIGrid3d g3d( 0, 1, 0, 1, 0, 1, 2, 6, 4, 6, comm3d);
        auto vec = dg::create::uniform_real_distribution( a, b, g3d);
        auto global = dg::create::uniform_real_distribution( a, b, g3d.global());
        auto local = dg::global2local( global, g3d);
        REQUIRE( vec.data().size() == local.data().size());
        for( unsigned u=0; u<vec.data().size(); u++)
        {
            INFO( "Local index "<<u);
            CHECK( vec.data()[u] == local.data()[u]);
        }
    }
}
