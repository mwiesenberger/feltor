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
    int status = 0;
    dg::exblas::reduce_mpi_cpu( 1, in.data(), out.data(), &status, comm2,
        comm_mod, comm_red);
    CHECK( dg::exblas::cpu::Round( out.data()) == size2);
    CHECK( status == 0);
    MPI_Comm_free( &comm2);
}

TEST_CASE( "Accumulate reduces the status with the superaccumulator")
{
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank( comm, &rank);
    MPI_Comm_size( comm, &size);
    MPI_Comm comm_mod, comm_red;
    dg::exblas::mpi_reduce_communicator( comm, &comm_mod, &comm_red);
    // Test the first and the last rank as well as all ranks
    auto failing = GENERATE( as<std::string>{}, "none", "first", "last", "all");
    int status = 0;
    if( failing == "all" || (failing == "first" && rank == 0)
        || (failing == "last" && rank == size-1))
        status = 1;
    // Every rank contributes 1 and the sum must not be affected by the status
    std::vector<int64_t> in( dg::exblas::BIN_COUNT, 0);
    std::vector<int64_t> out( in);
    dg::exblas::cpu::Accumulate( in.data(), 1.);
    dg::exblas::reduce_mpi_cpu( 1, in.data(), out.data(), &status, comm,
        comm_mod, comm_red);
    INFO( "Failing ranks: "<<failing);
    CHECK( dg::exblas::cpu::Round( out.data()) == size);
    CHECK( status == (failing == "none" ? 0 : 1));
}
