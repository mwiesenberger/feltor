#include <iostream>
#include <mpi.h>

#include "mpi_init.h"
#include "mpi_kron.h"
#include "exblas/mpi_accumulate.h"
#include "catch2/catch_test_macros.hpp"

// Count communicators that are freed and carry the count_free_keyval attribute
static int freed_count = 0;
static int count_free( MPI_Comm, int, void*, void*)
{
    freed_count++;
    return MPI_SUCCESS;
}

TEST_CASE("MPI Kron test")
{
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    SECTION( "Sub returns the cached communicator")
    {
        std::vector<int> np3d = {0,0,0};
        MPI_Dims_create( size, 3, &np3d[0]);
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, np3d,
            {false,false,false}, true);

        MPI_Comm comm_sub100 = dg::mpi_cart_sub( comm, {1,0,0});
        int ndims;
        dg::mpi_cartdim_get( comm_sub100, &ndims);
        CHECK( ndims == 1);
        int coords[1], np[1] = {0}, periods[1];
        dg::mpi_cart_get( comm_sub100, 1, np, periods, coords);
        CHECK( np[0] == np3d[0]);

        CHECK( dg::mpi_cart_sub( comm, {1,0,0}) == comm_sub100);
        // Subbing all dimensions returns the communicator itself
        CHECK( dg::mpi_cart_sub( comm, {1,1,1}) == comm);
        // A sub of a sub is the same as the sub of the root
        MPI_Comm comm_sub101 = dg::mpi_cart_sub( comm, {1,0,1});
        CHECK( dg::mpi_cart_sub( comm_sub101, {1,0}) == comm_sub100);
        CHECK( dg::mpi_cart_sub( comm_sub101, {0,1}) ==
            dg::mpi_cart_sub( comm, {0,0,1}));
        MPI_Comm_free( &comm);
    }
    SECTION( "Wrong number of dimensions throws")
    {
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0},
            {false,false,false}, true);
        CHECK_THROWS_AS( dg::mpi_cart_sub( comm, {1,0}), dg::Error);
        CHECK_THROWS_AS( dg::mpi_cart_sub( comm, {1,0,0,0}), dg::Error);
        CHECK_THROWS_AS( dg::mpi_cart_split_as<4>( comm), dg::Error);
        CHECK( dg::mpi_cart_split_as<2>( comm).size() == 2);
        MPI_Comm_free( &comm);
    }
    SECTION( "Check double sub does not segfault")
    {
        // This tests for issue in openmpi 4.0.0
        // https://github.com/open-mpi/ompi/issues/13081
        // Test passes if no segfault
        int rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);

        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {1,1});
        auto comms = dg::mpi_cart_split_as<2>(comm);
        MPI_Comm comm_join = dg::mpi_cart_kron(comms);
        MPI_Comm comm_mod, comm_red;
        dg::exblas::mpi_reduce_communicator( comm_join, &comm_mod, &comm_red);
        MPI_Comm comm2 = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {1,1,1});
        auto comms2 = dg::mpi_cart_split_as<2>(comm2);

        CHECK( comms.size() == 2);
        CHECK( comms2.size() == 2);
        MPI_Comm_free( &comm);
        MPI_Comm_free( &comm2);
    }

    SECTION( "Direct kronecker")
    {
        std::vector<int> np3d = {0,0,0};
        MPI_Dims_create( size, 3, &np3d[0]);
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, np3d, {0,0,0},
            true);
        MPI_Comm comm_sub100 = dg::mpi_cart_sub( comm, {1,0,0});
        MPI_Comm comm_sub001 = dg::mpi_cart_sub( comm, {0,0,1});

        MPI_Comm kron = dg::mpi_cart_kron( {comm_sub100, comm_sub001}); // 100 + 001
        int coords[2], np[2] = {0,0}, periods[2];
        dg::mpi_cart_get( kron, 2, np, periods, coords);
        CHECK( np[0] == np3d[0]);
        CHECK( np[1] == np3d[2]);
        CHECK( kron == dg::mpi_cart_sub( comm, {1,0,1}));
        MPI_Comm comm_sub010 = dg::mpi_cart_sub( comm, {0,1,0});
        MPI_Comm comm111 = MPI_COMM_NULL;
        CHECK_THROWS_AS(
            // Cannot create kronecker in reverse order
            comm111 = dg::mpi_cart_kron( {kron, comm_sub010}), // 101 + 010
            dg::Error
        );
        CHECK( comm111 == MPI_COMM_NULL);
        MPI_Comm doesNotWork = comm111;
        CHECK_THROWS_AS(
            doesNotWork =dg::mpi_cart_kron( {comm_sub010, comm_sub010}),
            dg::Error
        );
        CHECK( doesNotWork == MPI_COMM_NULL);
        std::array<MPI_Comm,3> axes = dg::mpi_cart_split_as<3>( comm);
        CHECK( axes[0] == comm_sub100);
        CHECK( axes[1] == comm_sub010);
        CHECK( axes[2] == comm_sub001);
        MPI_Comm joined = dg::mpi_cart_kron( axes);
        CHECK( joined == comm);

        // A test originating from geometry_elliptic_mpib
        auto taxe = dg::mpi_cart_kron( {axes[0], axes[1]});
        auto paxes = dg::mpi_cart_split_as<2>( taxe);
        auto kaxe = dg::mpi_cart_kron( paxes);
        CHECK( kaxe == taxe);
        MPI_Comm_free( &comm);
    }
    SECTION( "Kronecker of one or no communicator")
    {
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {0,0});
        MPI_Comm comm_sub10 = dg::mpi_cart_sub( comm, {1,0});
        CHECK( dg::mpi_cart_kron( {comm_sub10}) == comm_sub10);
        CHECK( dg::mpi_cart_kron( std::vector<MPI_Comm>{}) == MPI_COMM_NULL);
        MPI_Comm_free( &comm);
    }
    SECTION( "Kronecker of communicators with different roots throws")
    {
        MPI_Comm comm0 = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {0,0});
        MPI_Comm comm1 = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {0,0});
        MPI_Comm comm0_sub10 = dg::mpi_cart_sub( comm0, {1,0});
        MPI_Comm comm1_sub01 = dg::mpi_cart_sub( comm1, {0,1});
        CHECK_THROWS_AS( dg::mpi_cart_kron( {comm0_sub10, comm1_sub01}),
            dg::Error);
        MPI_Comm_free( &comm0);
        MPI_Comm_free( &comm1);
    }
    SECTION( "Kronecker of communicators unknown to dg throws")
    {
        MPI_Comm comm0 = dg::mpi_cart_create( MPI_COMM_WORLD, {0}, {0});
        MPI_Comm comm1 = dg::mpi_cart_create( MPI_COMM_WORLD, {0}, {0});
        CHECK_THROWS_AS( dg::mpi_cart_kron( {comm0, comm1}), dg::Error);
        MPI_Comm_free( &comm0);
        MPI_Comm_free( &comm1);
    }
    SECTION( "Freed communicators are not confused with new ones")
    {
        // MPI may reuse the handles of freed communicators
        int count_free_keyval;
        MPI_Comm_create_keyval( MPI_COMM_NULL_COPY_FN, count_free,
            &count_free_keyval, nullptr);
        for( unsigned i=0; i<4; i++)
        {
            MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {0,0});
            MPI_Comm comm_sub10 = dg::mpi_cart_sub( comm, {1,0});
            MPI_Comm comm_sub01 = dg::mpi_cart_sub( comm, {0,1});
            int ndims;
            dg::mpi_cartdim_get( comm_sub10, &ndims);
            CHECK( ndims == 1);
            dg::mpi_cartdim_get( comm_sub01, &ndims);
            CHECK( ndims == 1);
            CHECK( comm_sub10 != comm_sub01);
            CHECK( dg::mpi_cart_kron( {comm_sub10, comm_sub01}) == comm);
            MPI_Comm_set_attr( comm_sub10, count_free_keyval, nullptr);
            MPI_Comm_set_attr( comm_sub01, count_free_keyval, nullptr);
            // Freeing a cached communicator before the root must not lead to
            // a double free
            if( i%2 == 0)
                MPI_Comm_free( &comm_sub10);
            freed_count = 0;
            MPI_Comm_free( &comm);
            // The root frees the remaining cached communicators
            CHECK( freed_count == ( i%2 == 0 ? 1 : 2));
        }
        MPI_Comm_free_keyval( &count_free_keyval);
    }
    SECTION( "Documentation")
    {
        //! [split and join]
        // Create a Cartesian communicator
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {0,0,0});
        // Split into 3 axes
        std::array<MPI_Comm,3> axes3d = dg::mpi_cart_split_as<3>( comm);
        // Join the first and third axis
        MPI_Comm comm_101 = dg::mpi_cart_kron( {axes3d[0], axes3d[2]});
        // Split up again
        std::array<MPI_Comm,2> axes2d = dg::mpi_cart_split_as<2>( comm_101);
        CHECK( axes2d[0] == axes3d[0]);
        CHECK( axes2d[1] == axes3d[2]);
        //! [split and join]
        MPI_Comm_free( &comm);
    }
}

