/**
 *  @file mpi_accumulate.h
 *  @brief Primitives for an MPI Reduction
 *
 *  @authors
 *        Matthias Wiesenberger -- mattwi@fysik.dtu.dk
 */
#pragma once
#include <mpi.h>
#include <algorithm>
#include <array>
#include <vector>
#include "accumulate.h"

namespace dg
{
namespace exblas {

///@cond
namespace detail
{
// Free the cached communicators when the communicator they belong to is freed
inline int free_reduce_communicators( MPI_Comm, int, void* attribute, void*)
{
    auto* comms = static_cast<std::array<MPI_Comm, 2>*>( attribute);
    MPI_Comm_free( &(*comms)[0]);
    MPI_Comm_free( &(*comms)[1]);
    delete comms;
    return MPI_SUCCESS;
}
}//namespace detail
///@endcond

/**
 * @brief This function can be used to partition communicators for the \c
 * exblas::reduce_mpi_cpu function
 *
 * If the ranks in \c comm are aligned in rows of 128
 * \c comm_mod is the 1d communicator for the rows and \c comm_mod_reduce
 * is the communicator for the columns. The last row may contain fewer than 128
 * ranks, so the last columns of \c comm_mod_reduce may contain 1 rank fewer
 * than the others
 *
 * @param comm the input communicator (unmodified, may not be \c MPI_COMM_NULL)
 * @param comm_mod the line communicator in \c comm, consists of all rank/mod
 * ranks
 * @param comm_mod_reduce the column communicator in \c comm consists of all
 * rank%mod ranks
 * @attention The creation of new communicators involves communication between all
 * participation processes (comm in this case).
 * In order to avoid excessive creation of new MPI communicators
 * (there is a limit to how many a program can create), the function caches
 * the created communicators as an attribute of \c comm. If you
 * repeatedly call this function with the same \c comm only the first call will
 * actually create new communicators. They are freed when \c comm is freed and
 * are not copied by \c MPI_Comm_dup.
 */
inline void mpi_reduce_communicator(MPI_Comm comm,
    MPI_Comm* comm_mod, MPI_Comm* comm_mod_reduce){
    // We cache the communicators as an attribute of comm: unlike a map keyed
    // by the handle this is not fooled by a freed and reused handle and all
    // ranks agree on whether the communicators exist
    // The communicators are not copied by MPI_Comm_dup (MPI_COMM_NULL_COPY_FN):
    // a duplicate creates its own on first use, which keeps the separate
    // communication context of the duplicate and avoids collective calls in
    // the copy callback (that is also called by the non-blocking MPI_Comm_idup)
    static int keyval = [](){
        int key;
        MPI_Comm_create_keyval( MPI_COMM_NULL_COPY_FN,
            detail::free_reduce_communicators, &key, nullptr);
        return key;
    }();
    assert( comm != MPI_COMM_NULL);
    // The attribute is a pointer: MPI_Comm_set_attr takes the pointer itself,
    // but MPI_Comm_get_attr takes the address of a pointer (declared as void*
    // although it is used as void**), into which it writes the stored pointer
    std::array<MPI_Comm, 2>* comms;
    int found;
    MPI_Comm_get_attr( comm, keyval, &comms, &found);
    if( found)
    {
        *comm_mod = (*comms)[0];
        *comm_mod_reduce = (*comms)[1];
        return;
    }
    int rank, size;
    MPI_Comm_rank( comm, &rank);
    MPI_Comm_size( comm, &size);
    // For example
    // xxxxxxxxxxxxxxxx // 16 (comm_mod)
    // xxxxxxxxxxxxxxxx // 16
    // xxxxxx           //  6
    // 3333332222222222 // (comm_mod_reduce)
    int mod = 128;
    // 0-127, 128-255, ...
    MPI_Comm_split( comm, rank/mod, rank%mod, comm_mod); //collective call
    // Here we split rank%mod instead of returning MPI_COMM_NULL for != 0
    // https://github.com/open-mpi/ompi/issues/13081
    MPI_Comm_split( comm, rank%mod, rank, comm_mod_reduce);

    MPI_Comm_set_attr( comm, keyval,
        new std::array<MPI_Comm, 2>{*comm_mod, *comm_mod_reduce});
}

/*! @brief reduce a number of superaccumulators distributed among mpi processes

 * The words of a normalized superaccumulator are smaller than \f$ 2^{56}\f$
 * (\c exblas::DIGITS bits), so we cannot sum more than 128 accumulators in 64
 * bit integers before we need to normalize again. We therefore split the
 * reduction into several steps if more than 128 processes are involved.
 * This function normalizes, reduces, normalizes, reduces and
 * broadcasts the result to all participating processes.  As usual the resulting
 * superaccumulator is unnormalized.
 * @attention The second step sums one accumulator per row of 128 processes, so
 * the number of processes in \c comm must not exceed \f$ 128^2 = 16384\f$
 * @param num_superacc number of Superaccumulators each process holds
 * @param in unnormalized input superaccumulators ( must be of size
 * num_superacc*\c exblas::BIN_COUNT, allocated on the cpu) (read/write,
 * undefined on out)
 * @param out each process contains the result on output( must be of size
 * num_superacc*\c exblas::BIN_COUNT, allocated on the cpu) (write, may not
 * alias in)
 * @param status (read/write) on input the status of this process (0 for
 * success), on output 0 on every process if all processes had 0 and 1 on
 * every process otherwise
 * @param comm The complete MPI communicator
 * @param comm_mod This is the line communicator of up to 128 ranks
 * @param comm_mod_reduce This is the column communicator consisting of all
 * rank 0, 1, 2, ...,127 processes in all comm_mod
 * @sa \c exblas::mpi_reduce_communicator to generate the required communicators
*/
inline void reduce_mpi_cpu(  unsigned num_superacc, int64_t* in, int64_t* out,
int* status, MPI_Comm comm, MPI_Comm comm_mod, MPI_Comm comm_mod_reduce )
{
    for( unsigned i=0; i<num_superacc; i++)
    {
        int imin=exblas::IMIN, imax=exblas::IMAX;
        cpu::Normalize(&in[i*exblas::BIN_COUNT], imin, imax);
    }
    // The status travels in one additional word of the same message, which
    // saves a separate collective for every dot product. Only its sum matters
    // and it is never normalized
    const unsigned num_words = num_superacc*exblas::BIN_COUNT;
    std::vector<int64_t> send( in, in+num_words), receive( num_words+1);
    send.push_back( *status);
    int size_comm, size_comm_mod;
    MPI_Comm_size( comm, &size_comm);
    MPI_Comm_size( comm_mod, &size_comm_mod);
    if( size_comm == size_comm_mod) // size < 128
        MPI_Allreduce( send.data(), receive.data(), num_words+1, MPI_INT64_T, MPI_SUM, comm);
    else
    {
        MPI_Reduce( send.data(), receive.data(), num_words+1, MPI_INT64_T,
            MPI_SUM, 0, comm_mod);
        int rank;
        MPI_Comm_rank( comm_mod, &rank);
        if(rank == 0)
        {
            int size;
            MPI_Comm_size( comm_mod_reduce, &size);
            if( size > 1)
            {
                for( unsigned i=0; i<num_superacc; i++)
                {
                    int imin=exblas::IMIN, imax=exblas::IMAX;
                    cpu::Normalize(&receive[i*exblas::BIN_COUNT], imin, imax);
                }
                send = receive;
                MPI_Reduce( send.data(), receive.data(), num_words+1, MPI_INT64_T,
                    MPI_SUM, 0, comm_mod_reduce);
            }
        }
        MPI_Bcast( receive.data(), num_words+1, MPI_INT64_T, 0, comm);
    }
    std::copy_n( receive.data(), num_words, out);
    *status = receive[num_words] != 0 ? 1 : 0;
}

}//namespace exblas
} //namespace dg
