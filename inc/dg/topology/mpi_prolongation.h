#pragma once

#include "dg/backend/mpi_init.h"
#include "mpi_projection.h"
#include "prolongation.h"

/*! @file
  @brief Useful MPI typedefs and overloads of interpolation and projection
  */

namespace dg
{
namespace create
{
///@copydoc dg::create::prolongation(const aRealTopology<real_type,Nd>&,std::array<unsigned,Md>)
///@note In MPI there is no communication involved
template<class MPITopology, size_t Md, typename = std::enable_if_t< dg::is_vector_v<
    typename MPITopology::host_vector, MPIVectorTag>>>
dg::MIHMatrix_t<typename MPITopology::value_type> prolongation(
    const MPITopology& g_new, std::array<unsigned,Md> axes)
{
    using real_type = typename MPITopology::value_type;
    dg::IHMatrix_t<real_type> mat = dg::create::prolongation(
        g_new.local(), axes); // local rows, local cols
    return dg::MIHMatrix_t<real_type>( mat);
}

///@cond
namespace detail
{
// Return a gather that copies a vector from the process distribution of
// g_from into the process distribution of g_to (same global grid). Element i
// of the local result is element bufferIdx[i] of the gather's buffer.
//
// For every element of its local part of g_to, each process finds the
// process that holds this element in g_from. That process is given as a rank
// in the caller's g_from.communicator(), and the function translates it into
// a rank in comm. So comm must contain all processes of
// g_from.communicator(), all processes of comm must call this function, and
// comm must outlive the gather.
//
// The two calls in mpi_reduction, with comm = g.communicator() (P processes):
// - Forward, g_from = g and g_to = g2: both communicators have all P
//   processes, arranged in different process grids.
// - Back, g_from = g2.grid(remains) and g_to = g.grid(remains):
//   g_from.communicator() has all P processes (g2 has one process along the
//   averaged axes). g_to.communicator() has only the P/P_A processes that
//   share the caller's position along the averaged axes; each of these P_A
//   communicators receives the whole reduced vector. comm cannot be one of
//   the two grid communicators: g_to's lacks most owners and g_from's is
//   freed after the construction.
template<class MPITopologyTo, class MPITopologyFrom>
MPIGather<thrust::host_vector> make_redistribution( const MPITopologyTo& g_to,
    const MPITopologyFrom& g_from, MPI_Comm comm,
    thrust::host_vector<int>& bufferIdx)
{
    static_assert( MPITopologyTo::ndim() == MPITopologyFrom::ndim());
    assert( g_to.global().get_shape() == g_from.global().get_shape());
    // The ranks of g_from.communicator() as ranks in comm
    int from_size;
    MPI_Comm_size( g_from.communicator(), &from_size);
    std::vector<int> from_ranks( from_size), comm_ranks( from_size);
    for( int r=0; r<from_size; r++)
        from_ranks[r] = r;
    MPI_Group from_group, comm_group;
    MPI_Comm_group( g_from.communicator(), &from_group);
    MPI_Comm_group( comm, &comm_group);
    MPI_Group_translate_ranks( from_group, from_size, from_ranks.data(),
        comm_group, comm_ranks.data());
    MPI_Group_free( &from_group);
    MPI_Group_free( &comm_group);

    int rank_to;
    MPI_Comm_rank( g_to.communicator(), &rank_to);
    thrust::host_vector<std::array<int,2>> gather_map( g_to.local_size());
    for( int i=0; i<(int)gather_map.size(); i++)
    {
        int gIdx = 0, lIdx = 0, pid = 0;
        g_to.local2globalIdx( i, rank_to, gIdx);
        g_from.global2localIdx( gIdx, lIdx, pid);
        assert( comm_ranks[pid] != MPI_UNDEFINED);
        gather_map[i] = { comm_ranks[pid], lIdx};
    }
    return MPIGather<thrust::host_vector>( gather_map, bufferIdx, comm);
}

// The grid of g along the given axes
template<class MPITopology, size_t Md>
RealMPIGrid<typename MPITopology::value_type, Md> sub_grid(
    const MPITopology& g, std::array<unsigned,Md> axes)
{
    std::array<RealMPIGrid<typename MPITopology::value_type,1>,Md> gs;
    for( unsigned u=0; u<Md; u++)
        gs[u] = g.grid( axes[u]);
    return RealMPIGrid<typename MPITopology::value_type, Md>( gs);
}

// Same global grid as g, with one process along the given axes. The
// processes along the axes are moved to the remaining axes, outermost
// (largest axis number, slowest varying in memory) first, such that every
// process keeps at least one cell. Splitting the outermost axis keeps the
// messages of the gathers contiguous and the received elements in local
// order. The caller must free the communicator of the returned grid (this
// frees its sub-communicators).
template<class MPITopology, size_t Md>
RealMPIGrid<typename MPITopology::value_type, MPITopology::ndim()>
    redistributed_grid( const MPITopology& g, std::array<unsigned,Md> axes)
{
    constexpr size_t Nd = MPITopology::ndim();
    std::vector<int> dims( Nd), periods( Nd), coords( Nd);
    dg::mpi_cart_get( g.communicator(), Nd, &dims[0], &periods[0], &coords[0]);
    int factor = 1;
    for( unsigned u=0; u<Md; u++)
    {
        factor *= dims[axes[u]];
        dims[axes[u]] = 1;
    }
    std::array<unsigned,Nd-Md> remains = detail::complement( g, axes);
    // Try all divisors of the factor, the largest first, on the outermost
    // remaining axis, then recurse to the next inner one
    auto distribute = [&]( auto& self, int k, int P) -> bool
    {
        if( P == 1)
            return true;
        if( k < 0)
            return false;
        unsigned u = remains[k];
        for( int f = P; f >= 1; f--)
        {
            if( P%f != 0 || (unsigned)(dims[u]*f) > g.N(u))
                continue;
            dims[u] *= f;
            if( self( self, k-1, P/f))
                return true;
            dims[u] /= f;
        }
        return false;
    };
    if( !distribute( distribute, int(Nd-Md)-1, factor))
    {
        int size;
        MPI_Comm_size( g.communicator(), &size);
        dg::Message msg(_ping_);
        msg << "Cannot distribute "<<size<<" processes over the remaining axes"
            <<" with cell numbers";
        for( unsigned u=0; u<Nd-Md; u++)
            msg << " "<<g.N(remains[u]);
        msg <<" such that every process holds at least one cell";
        throw dg::Error( msg);
    }
    MPI_Comm comm = dg::mpi_cart_create( g.communicator(), dims, periods);
    const auto& gg = g.global();
    return { gg.get_p(), gg.get_q(), gg.get_n(), gg.get_N(), gg.get_bc(),
        dg::mpi_cart_split_as<Nd>( comm)};
}

// Rows of the reduction (or projection) matrix of the local grid of g, with
// the binary exact weights of the MPI grid (local rows, local cols)
template<class MPITopology, size_t Md>
dg::IHMatrix_t<typename MPITopology::value_type> local_reduction(
    std::array<unsigned,Md> axes, const MPITopology& g, bool project)
{
    auto mat = dg::create::prolongation( g.local(), axes).transpose();
    if( project)
    {
        auto g_new = sub_grid( g, detail::complement( g, axes));
        detail::scale_projection( mat, dg::create::inv_weights( g_new).data(),
            dg::create::weights( g).data());
    }
    return mat;
}

// The reduction (or projection) as an MPI matrix that computes every row on
// one process in the order of the shared memory matrix
template<class MPITopology, size_t Md>
dg::MIHMatrix_t<typename MPITopology::value_type> mpi_reduction(
    std::array<unsigned,Md> axes, const MPITopology& g, bool project)
{
    using real_type = typename MPITopology::value_type;
    int factor = 1;
    for( unsigned u=0; u<Md; u++)
    {
        int size;
        MPI_Comm_size( g.comm( axes[u]), &size);
        factor *= size;
    }
    // The averaged axes are not distributed: all rows are local
    if( factor == 1)
        return dg::MIHMatrix_t<real_type>( local_reduction( axes, g, project));

    auto g2 = redistributed_grid( g, axes);
    auto remains = detail::complement( g, axes);
    thrust::host_vector<int> forward_idx, back_idx;
    auto forward = make_redistribution( g2, g, g.communicator(), forward_idx);
    auto back = make_redistribution( sub_grid( g, remains),
        sub_grid( g2, remains), g.communicator(), back_idx);
    // The local matrix acts on the buffer of forward
    auto mat = local_reduction( axes, g2, project);
    thrust::host_vector<int> cols = mat.column_indices();
    for( auto& c : cols)
        c = forward_idx[c];
    MPI_Comm g2_comm = g2.communicator();
    MPI_Comm_free( &g2_comm);
    return { dg::IHMatrix_t<real_type>( mat.num_rows(), forward.buffer_size(),
        mat.row_offsets(), cols, mat.values()), forward, back, back_idx};
}
}//namespace detail
///@endcond

/*! @class hide_mpi_reduction
 * @note In MPI, if an axis in \c axes is distributed among more than one
 * process, the matrix first gathers the input vector into a distribution with
 * one process along \c axes; the processes along \c axes are moved to the
 * remaining axes, outermost first. Each process then computes its rows in the
 * order of the shared memory matrix, and the result is gathered back to all
 * processes along \c axes. The result is therefore binary reproducible: it
 * is the same for every process distribution and equal to the shared memory
 * result on the same backend.
 * @throw dg::Error if the remaining axes have too few cells to give every
 * process at least one cell
 */

///@copydoc dg::create::reduction(std::array<unsigned,Md>,const aRealTopology<real_type,Nd>&)
///@copydoc hide_mpi_reduction
template<class MPITopology, size_t Md, typename = std::enable_if_t< dg::is_vector_v<
    typename MPITopology::host_vector, MPIVectorTag>>>
dg::MIHMatrix_t<typename MPITopology::value_type> reduction(
    std::array<unsigned,Md> axes, const MPITopology& g_old)
{
    return detail::mpi_reduction( axes, g_old, false);
}

///@copydoc dg::create::projection(std::array<unsigned,Md>,const aRealTopology<real_type,Nd>&)
///@copydoc hide_mpi_reduction
template<class MPITopology, size_t Md, typename = std::enable_if_t< dg::is_vector_v<
    typename MPITopology::host_vector, MPIVectorTag>>>
dg::MIHMatrix_t<typename MPITopology::value_type> projection(
    std::array<unsigned,Md> axes, const MPITopology& g_old)
{
    return detail::mpi_reduction( axes, g_old, true);
}
}//namespace create
}//namespace dg
