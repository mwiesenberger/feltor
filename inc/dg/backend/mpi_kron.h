#pragma once
#include <cassert>
#include <array>
#include <vector>
#include <algorithm> // std::copy_n, std::all_of
#include <map>
#include <memory>
#include "mpi_cart.h"
#include "exceptions.h"


namespace dg
{
/*! @class hide_mpi_cart_rationale
 * @note The reason the dg library provides \c dg::mpi_cart_sub and \c
 * dg::mpi_cart_kron is that unfortunately the MPI standard does not provide a
 * way to form the Kronecker product of Cartesian communicators without
 * manually tracking their parent Cartesian communicators.  However, this is
 * needed in the \c dg::blas1::kronecker and \c dg::kronecker functions.
 */
/*! @class hide_mpi_cart_ownership
 * @note The dg library caches the communicators that \c dg::mpi_cart_sub
 * creates. Repeat calls that ask for the same dimensions of the same root
 * communicator, i.e. the Cartesian communicator from which they derive and
 * that was created by \c dg::mpi_cart_create, return the same communicator.
 * All derived communicators are freed when the root is freed and must not be
 * freed by the user. The user must only free the root communicator.
 */
///@cond
namespace detail{
// The communicators that derive from one root Cartesian communicator through
// dg::mpi_cart_sub, keyed by their remain_dims relative to the root.
// The root itself is the entry with all remain_dims true
struct MPICartFamily
{
    std::map<std::vector<int>, MPI_Comm> members;
};
// The attribute of every member of a family
struct MPICartNode
{
    std::shared_ptr<MPICartFamily> family;
    std::vector<int> remain_dims; // relative to the root
};

// Called by MPI when a member is freed; the root frees all other members
inline int mpi_cart_delete_node( MPI_Comm comm, int, void* attribute, void*)
{
    auto* node = static_cast<MPICartNode*>( attribute);
    auto& members = node->family->members;
    const std::vector<int>& re = node->remain_dims;
    if( std::all_of( re.begin(), re.end(), [](int r){ return r;}))
    {
        // Freeing a member calls this function, which erases the member from
        // members, so we cannot iterate over members and use a copy instead
        auto children = members;
        for( auto& child : children)
            if( child.second != comm)
                MPI_Comm_free( &child.second);
    }
    else
        members.erase( re);
    delete node;
    return MPI_SUCCESS;
}

// The attribute is not copied by MPI_Comm_dup (MPI_COMM_NULL_COPY_FN): a
// duplicate becomes the root of its own family on first use
inline int mpi_cart_keyval()
{
    static int keyval = [](){
        int key;
        MPI_Comm_create_keyval( MPI_COMM_NULL_COPY_FN, mpi_cart_delete_node,
            &key, nullptr);
        return key;
    }();
    return keyval;
}

// Return the node of comm or nullptr if comm is not a member of a family
inline const MPICartNode* mpi_cart_find( MPI_Comm comm)
{
    // MPI_Comm_get_attr takes the address of a pointer (declared as void*
    // although it is used as void**), into which it writes the stored pointer
    MPICartNode* node;
    int found;
    MPI_Comm_get_attr( comm, mpi_cart_keyval(), &node, &found);
    return found ? node : nullptr;
}

// Add comm to family
inline void mpi_cart_add( MPI_Comm comm, std::shared_ptr<MPICartFamily>
    family, const std::vector<int>& remain_dims)
{
    family->members[remain_dims] = comm;
    MPI_Comm_set_attr( comm, mpi_cart_keyval(),
        new MPICartNode{ family, remain_dims});
}

// Return the node of comm, make comm the root of a new family if it is not a
// member of a family yet
inline const MPICartNode& mpi_cart_node( MPI_Comm comm)
{
    const MPICartNode* node = mpi_cart_find( comm);
    if( node == nullptr)
    {
        int ndims;
        dg::mpi_cartdim_get( comm, &ndims);
        mpi_cart_add( comm, std::make_shared<MPICartFamily>(),
            std::vector<int>( ndims, true));
        node = mpi_cart_find( comm);
    }
    return *node;
}
}
///@endcond
///@addtogroup mpi_utility
///@{

/*! @brief Call \c MPI_Cart_sub or return the communicator of a previous call
 *
 * @param comm communicator with Cartesian structure (handle) (parameter used
 * in \c MPI_Cart_sub)
 * @param remain_dims the i-th entry of \c remain_dims specifies whether the
 * (N-1-i)-th dimension, which in dg ordering is the i-th dimension, is kept
 * in the subgrid (true) or is dropped (false), must
 * have \c ndims entries, else \c dg::Error is thrown. (parameter used in
 * \c MPI_Cart_sub)
 * @note Notice the transposition in \c remain_dims ! In the dg-library we
 * always transpose Cartesian communicators to effectively make the process
 * ranking column-major instead of row-major.
 * @return communicator containing the subgrid that includes the calling
 * process (handle) (parameter used in \c MPI_Cart_sub). If a previous call
 * already created the communicator with the same dimensions of the same root,
 * that communicator is returned.
 * @copydoc hide_mpi_cart_ownership
 * @copydoc hide_mpi_cart_rationale
 */
inline MPI_Comm mpi_cart_sub( MPI_Comm comm, std::vector<int> remain_dims)
{
    int ndims;
    dg::mpi_cartdim_get( comm, &ndims);
    if( (unsigned) ndims != remain_dims.size())
        throw Error(Message(_ping_)<<"remain_dims has "<<remain_dims.size()
                <<" entries but comm has "<<ndims<<" dimensions");

    const detail::MPICartNode& node = detail::mpi_cart_node( comm);
    // The true entries of node.remain_dims are the dimensions of comm
    // root_dims may be larger than remain_dims because comm may have a
    // parent but exactly remain_dims.size() entries are true
    std::vector<int> root_dims = node.remain_dims;
    unsigned counter = 0;
    for( auto& r : root_dims)
        if( r)
            r = (remain_dims[counter++] != 0);
    assert( counter == remain_dims.size());
    auto& members = node.family->members;
    auto it = members.find( root_dims);
    if( it != members.end())
        return it->second;
    MPI_Comm newcomm;
    int err = dg::mpi_cart_sub( comm, remain_dims.data(), &newcomm);
    if( err != MPI_SUCCESS)
        throw Error(Message(_ping_)<<
                "Cannot create Cartesian sub comm from given communicator");
    detail::mpi_cart_add( newcomm, node.family, root_dims);
    return newcomm;
}

/*! @brief Form a (transposed) Kronecker product among Cartesian communicators
 *
 * Find communicator as the one that hypothetically generated all
 * input comms through <tt>dg::mpi_cart_sub( return_comm, remain_dims[u],
 * comms[u]);</tt> for all <tt>u < comms.size();</tt>
 *
 * Unless \c comms is empty or contains only 1 communicator all input comms
 * must have been returned by \c dg::mpi_cart_sub (or \c dg::mpi_cart_split,
 * \c dg::mpi_cart_kron) from the same root Cartesian communicator, else
 * \c dg::Error is thrown.
 * Furthermore the comms must be mutually orthogonal i.e. any \c true entry in
 * \c remain_dims can exist in only exactly one comm.  The resulting \c
 * remain_dims of the output is then the union of all \c remain_dims of the
 * inputs.
 *
 * For example
 * @snippet{trimleft} mpi_kron_mpit.cpp split and join
 *
 * @attention The order of communicators matters. The function will always
 * transpose communicators, such that the first dimension in \c comms is the
 * last dimension in the return Cartesian communicator
 * @param comms input communicators
 * @return Kronecker product of communicators. If \c comms is empty \c
 * MPI_COMM_NULL is returned, if it contains only one communicator that
 * communicator is returned
 * @copydoc hide_mpi_cart_ownership
 * @copydoc hide_mpi_cart_rationale
 */
inline MPI_Comm mpi_cart_kron( std::vector<MPI_Comm> comms)
{
    // This non-template interface must exist so compiler can deduce call
    // mpi_cart_kron( {comm0, comm1});
    if ( comms.empty())
        return MPI_COMM_NULL;
    if( comms.size() == 1)
        return comms[0];

    std::vector<const detail::MPICartNode*> nodes( comms.size());
    for( unsigned u=0; u<comms.size(); u++)
    {
        nodes[u] = detail::mpi_cart_find( comms[u]);
        if( nodes[u] == nullptr)
            throw Error(Message(_ping_)<<
                    "In mpi_cart_kron comm number "<<u
                    <<" was not created by dg::mpi_cart_sub");
        // this just compares the two std::shared_ptr -> very cheap: 
        if( nodes[u]->family != nodes[0]->family)
            throw Error(Message(_ping_)<<
                    "In mpi_cart_kron all comms must have same root comm. "
                    <<"Offending comm number "<<u);
    }
    size_t ndims = nodes[0]->remain_dims.size();
    std::vector<int> remain_dims( ndims, false) ;
    unsigned current_free_k=0;
    for( unsigned u=0; u<comms.size(); u++)
    {
        for( unsigned k=0; k<ndims; k++)
        if( nodes[u]->remain_dims[k])
        {
            if( remain_dims[k])
                throw Error(Message(_ping_)<<
                    "Cannot form kronecker product with given communicators as remaining dims must be orthogonal ");
            if( k < current_free_k)
                throw Error(Message(_ping_)<<
                    "Cannot form kronecker product with communicators in reverse order compared to original ");
            remain_dims[k] = nodes[u]->remain_dims[k];
            current_free_k = k+1;
        }
    }
    MPI_Comm root = nodes[0]->family->members.at( std::vector<int>( ndims,
        true));
    return dg::mpi_cart_sub( root, remain_dims);
}

/*!
 * @brief Convenience shortcut for <tt> return mpi_cart_kron(
 * std::vector<MPI_Comm>(comms.begin(), comms.end()));</tt>
 */
template<class Vector>
MPI_Comm mpi_cart_kron( const Vector& comms)
{
    return mpi_cart_kron( std::vector<MPI_Comm>(comms.begin(), comms.end()));
}


/*! @brief Split a Cartesian communicator along each dimensions
 *
 * using repeated calls to \c dg::mpi_cart_sub
 * @param comm input Cartesian communicator
 * @return Array of 1-dimensional Cartesian communicators. The i-th
 * communicator corresponds to the (N-1-i)-th axis, which in dg ordering is
 * the i-th axis, in \c comm
 * @attention Notice the transpose. The first dimension in \c comm is the last
 * dimension in the return vector.
 * @copydoc hide_mpi_cart_ownership
 */
inline std::vector<MPI_Comm> mpi_cart_split( MPI_Comm comm)
{
    int ndims;
    dg::mpi_cartdim_get( comm, &ndims);

    std::vector<MPI_Comm> comms(ndims);
    for( int u=0; u<ndims; u++)
    {
        std::vector<int> remain_dims(ndims, 0);
        remain_dims[u] = 1;
        comms[u] = mpi_cart_sub( comm, remain_dims);
    }
    return comms;
}
/*!
 * @brief Same as \c mpi_cart_split but different return type
 *
 * @snippet{trimleft} mpi_kron_mpit.cpp split and join
 * @tparam Nd Number of dimensions to copy from \c mpi_cart_split
 * @param comm input Cartesian communicator ( <tt>Nd <= comm.ndims</tt>, else
 * \c dg::Error is thrown)
 * @return Array of 1-dimensional Cartesian communicators
 * @attention Notice the transpose. The first dimension in \c comm is the last
 * dimension in the return vector.
 * @copydoc hide_mpi_cart_ownership
 */
template<size_t Nd>
std::array<MPI_Comm, Nd> mpi_cart_split_as( MPI_Comm comm)
{
    int ndims;
    dg::mpi_cartdim_get( comm, &ndims);
    if( Nd > (unsigned)ndims)
        throw Error(Message(_ping_)<<"Cannot copy "<<Nd<<" axes from comm with "
                <<ndims<<" dimensions");
    auto split = mpi_cart_split( comm);
    std::array<MPI_Comm, Nd> arr;
    std::copy_n( split.begin(), Nd, arr.begin());
    return arr;
}



///@}


}
