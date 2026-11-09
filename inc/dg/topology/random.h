#pragma once

#include <thrust/host_vector.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/uniform_real_distribution.h>
#include "grid.h"
#ifdef MPI_VERSION
#include "dg/backend/mpi_vector.h"
#include "mpi_grid.h"
#endif //MPI_VERSION

/*! @file
  @brief Random vectors that do not depend on the parallelization
  */
namespace dg
{
namespace create
{
///@cond
namespace detail
{
template<class real_type, size_t Nd>
thrust::host_vector<real_type> do_uniform_real_distribution( real_type a,
    real_type b, const aRealTopology<real_type,Nd>& g)
{
    thrust::minstd_rand rng;
    thrust::uniform_real_distribution<real_type> dist( a, b);
    thrust::host_vector<real_type> vec( g.size());
    for( unsigned i=0; i<vec.size(); i++)
        vec[i] = dist( rng);
    return vec;
}
#ifdef MPI_VERSION
template<class real_type, size_t Nd>
MPI_Vector<thrust::host_vector<real_type>> do_uniform_real_distribution(
    real_type a, real_type b, const aRealMPITopology<real_type,Nd>& g)
{
    thrust::uniform_real_distribution<real_type> dist( a, b);
    thrust::host_vector<real_type> vec( g.local().size());
    int rank;
    MPI_Comm_rank( g.communicator(), &rank);
    for( unsigned idx=0; idx<vec.size(); idx++)
    {
        int gIdx = 0;
        if( !g.local2globalIdx( idx, rank, gIdx))
            throw Error( Message(_ping_)<<"Local index "<<idx<<" on rank "<<rank
                <<" has no global int index! Global grid size is "<<g.size());
        // The distribution draws exactly one number from the engine
        thrust::minstd_rand rng;
        rng.discard( gIdx);
        vec[idx] = dist( rng);
    }
    return MPI_Vector<thrust::host_vector<real_type>>( vec, g.communicator());
}
#endif //MPI_VERSION
}//namespace detail
///@endcond

///@addtogroup evaluation
///@{

/*! @brief Uniformly distributed random numbers on a grid that do not depend
 * on the parallelization
 *
 * The element with (global) index \c i of the returned vector is the
 * <tt>i</tt>-th number drawn with
 * <tt>thrust::uniform_real_distribution<real_type>( a, b)</tt> from a
 * default constructed <tt>thrust::minstd_rand</tt> engine. The result is
 * thus bitwise identical for a shared memory grid and an MPI grid with the
 * same global size, independent of the number of processes.
 *
 * Use this function instead of \c dg::RandomNumbers, which yields results
 * that depend on the number of threads (and processes).
 * @code{.cpp}
 * dg::x::DVec rnd = dg::construct<dg::x::DVec>(
 *     dg::create::uniform_real_distribution( 0, 1, grid));
 * @endcode
 * @tparam Topology A shared memory (derived from \c dg::aRealTopology) or MPI
 * (derived from \c dg::aRealMPITopology) grid
 * @param a lower bound of the interval
 * @param b upper bound of the interval
 * @param g The grid
 * @return Random numbers in <tt>[a,b)</tt> as a host vector (shared memory)
 * or an MPI vector of host vectors (MPI). Its value type is the
 * same as the grid value type
 * @note The MPI version computes each element separately with
 * <tt>thrust::minstd_rand::discard</tt>, which is \f$O(\log i)\f$, so the
 * cost is \f$ O( N \log N)\f$ for \f$N\f$ local points
 */
template<class Topology>
auto uniform_real_distribution( typename Topology::value_type a,
    typename Topology::value_type b, const Topology& g)
{
    return detail::do_uniform_real_distribution( a, b, g);
}
///@}

}//namespace create
}//namespace dg
