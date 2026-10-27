#include <iostream>
#include "exblas/mpi_accumulate.h"



#include "catch2/catch_all.hpp"

// TODO Find a test to test accumulation for >128 threads
TEST_CASE( "Accumulate reuses comm")
{
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank( comm, &rank);

    MPI_Comm comm_mod, comm_red;
    dg::exblas::mpi_reduce_communicator( comm, &comm_mod, &comm_red);
    // No new communicators should be generated
    MPI_Comm comm_mod2, comm_red2;
    dg::exblas::mpi_reduce_communicator( comm, &comm_mod2, &comm_red2);
    dg::exblas::mpi_reduce_communicator( comm, &comm_mod2, &comm_red2);
    dg::exblas::mpi_reduce_communicator( comm, &comm_mod2, &comm_red2);
    CHECK( comm_mod == comm_mod2);
    CHECK( comm_red == comm_red2);

}

TEST_CASE( "Accumulate does not reuse comm of a freed comm")
{
    // MPI may reuse the handle of a freed communicator for a new one
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm comm1;
    MPI_Comm_dup( MPI_COMM_WORLD, &comm1);
    MPI_Comm comm_mod, comm_red;
    dg::exblas::mpi_reduce_communicator( comm1, &comm_mod, &comm_red);
    MPI_Comm_free( &comm1);

    MPI_Comm comm2;
    MPI_Comm_split( MPI_COMM_WORLD, rank%2, rank, &comm2);
    dg::exblas::mpi_reduce_communicator( comm2, &comm_mod, &comm_red);
    int size2, size_mod;
    MPI_Comm_size( comm2, &size2);
    MPI_Comm_size( comm_mod, &size_mod);
    CHECK( size_mod == size2); // fewer than 128 ranks

    // Every rank contributes 1
    std::vector<int64_t> in( dg::exblas::BIN_COUNT, 0), out( in);
    dg::exblas::cpu::Accumulate( in.data(), 1.);
    dg::exblas::reduce_mpi_cpu( 1, in.data(), out.data(), comm2, comm_mod,
        comm_red);
    CHECK( dg::exblas::cpu::Round( out.data()) == size2);
    MPI_Comm_free( &comm2);
}
