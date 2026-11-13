#ifndef _DG_VECTOR_TRAITS_
#define _DG_VECTOR_TRAITS_

#include <type_traits>
#include <thrust/memory.h>
#include "scalar_categories.h"
#include "vector_categories.h"
#include "matrix_categories.h"
#include "execution_policy.h"

// Q: Make inclusion of tensor_traits easier if all you want is to use the full traits system in feltor?
// A: Some subtle issues:
// - Do you want the thrust header? (And thus depend on thrust library?)
// - Do you want all matrix traits as well?
// - If so, do you want the sparsematrix.h which incurs -lcusparse dependency
// In fact, it is possible to include tensor_traits.h as is and include the specialisation later when it is needed

namespace dg{
/*!@addtogroup traits
 *
 * @anchor dg_tensor_traits_table
 * The \c dg::TensorTraits specialisations of the dg library define
 *
 * | Type | `value_type` | `tensor_category` | `execution_policy` |
 * |------|--------------|-------------------|--------------------|
 * | any other type (primary template) | `void` | `dg::NotATensorTag` | `dg::NoPolicyTag` |
 * | `T` floating point (e.g. `float`, `double`) | `T` | `dg::FloatingPointTag` | `dg::AnyPolicyTag` |
 * | `T` integral (e.g. `int`, `bool`, `char`) | `T` | `dg::IntegralTag` | `dg::AnyPolicyTag` |
 * | `std::complex<T>` | `std::complex<T>` | `dg::ComplexTag` | `dg::AnyPolicyTag` |
 * | `thrust::complex<T>` | `thrust::complex<T>` | `dg::ComplexTag` | `dg::AnyPolicyTag` |
 * | `thrust::host_vector<T>` | `T` | `dg::ThrustVectorTag` | `dg::SerialTag` |
 * | `thrust::device_vector<T>` | `T` | `dg::ThrustVectorTag` | `dg::CudaTag`, `dg::OmpTag` or `dg::SerialTag` depending on `THRUST_DEVICE_SYSTEM` |
 * | `std::vector<T>`, `T` scalar | `T` | `dg::ThrustVectorTag` | `dg::SerialTag` |
 * | `std::vector<T>`, `T` not scalar | `dg::get_value_type<T>` | `dg::RecursiveVectorTag` | `dg::get_execution_policy<T>` |
 * | `std::array<T,N>`, `T` scalar | `T` | `dg::ArrayScalarTag` | `dg::SerialTag` |
 * | `std::array<T,N>`, `T` not scalar | `dg::get_value_type<T>` | `dg::ArrayVectorTag` | `dg::get_execution_policy<T>` |
 * | `std::map<Key,T>` | `dg::get_value_type<T>` | `dg::StdMapTag` | `dg::get_execution_policy<T>` |
 * | `dg::View<V>` | `dg::get_value_type<V>` | `dg::ThrustVectorTag` | `dg::get_execution_policy<V>` |
 * | `dg::MPI_Vector<V>` | `dg::get_value_type<V>` | `dg::MPIVectorTag` | `dg::get_execution_policy<V>` |
 * | `dg::EllSparseBlockMat<T,V>`, `dg::CooSparseBlockMat<T,V>` | `T` | `dg::SparseBlockMatrixTag` | `dg::get_execution_policy<V<T>>` |
 * | `dg::SparseMatrix<I,T,V>` | `T` | `dg::SparseMatrixTag` | `dg::SparseMatrix<I,T,V>::policy` |
 * | `dg::DenseMatrix<V>` | `dg::get_value_type<V>` | `dg::DenseMatrixTag` | `dg::get_execution_policy<V>` |
 * | `dg::MPISparseBlockMat<V,LI,LO>`, `dg::MPIDistMat<V,LI,LO>` | `dg::get_value_type<LI>` | `dg::MPIMatrixTag` | none |
 * | `dg::SquareMatrix<T>` | `T` | `dg::SelfMadeMatrixTag` | `dg::SerialTag` |
 * | `dg::Composite<M>` | `dg::get_value_type<M>` | `dg::SelfMadeMatrixTag` | none |
 * | `dg::Elliptic1d`, `dg::Elliptic2d`, `dg::Elliptic3d`, `dg::RefinedElliptic`, `dg::GeneralHelmholtz`, `dg::Helmholtz2`, `dg::ChebyshevPreconditioner`, `dg::ModifiedChebyshevPreconditioner`, `dg::LeastSquaresPreconditioner`, `dg::MultiMatrix` with container `V` | `dg::get_value_type<V>` | `dg::SelfMadeMatrixTag` | none |
 *
 * "none" means the specialisation does not define the member.
 */
///@addtogroup dispatch
///@{


/*! @brief The vector traits

Specialize this struct if you want to enable your own vector/container class for the use in blas1 functions.
The contained types are
- <tt> value_type </tt> the elementary data type of the contained data
- <tt> tensor_category </tt> the layout of the data (derives from \c AnyMatrixTag)
- <tt> execution_policy </tt> for \c SharedVectorTag the execution policy
    (derives from \c AnyPolicyTag)
    indicates the type of hardware memory is physically
allocated on in a vector class and therefore indicates the
possible parallelization and optimization strategies.

Matrices may omit \c execution_policy (the MPI matrices and most self-made
matrices do); \c dg::get_execution_policy is then ill-formed for them.
All specialisations are listed in the \ref dg_tensor_traits_table "TensorTraits table".
\see \ref dispatch
*/
template< class Vector, class Enable=void>
struct TensorTraits
{
    using value_type = void; //!< see the class description
    using tensor_category = NotATensorTag; //!< see the class description
    using execution_policy = NoPolicyTag; //!< see the class description
};

/// Shortcut for <tt>dg::TensorTraits<std::decay_t<Vector>>::value_type</tt>
template<class Vector>
using get_value_type = typename TensorTraits<std::decay_t<Vector>>::value_type;
/// Shortcut for <tt>dg::TensorTraits<std::decay_t<Vector>>::tensor_category</tt>
template<class Vector>
using get_tensor_category = typename TensorTraits< std::decay_t<Vector>>::tensor_category;
/// Shortcut for <tt>dg::TensorTraits<std::decay_t<Vector>>::execution_policy</tt>
template<class Vector>
using get_execution_policy = typename TensorTraits<std::decay_t<Vector>>::execution_policy;

///@}

///@cond
////////////get element, pointer and data
template<class T> //T = SharedVector
using get_pointer_type = std::conditional_t< std::is_const< std::remove_reference_t<T> >::value,
    const get_value_type<T>*, get_value_type<T>* >;

template<class T> //T = RecursiveVector
using get_element_type = std::conditional_t< std::is_const< std::remove_reference_t<T> >::value,
    const typename std::decay_t<T>::value_type&, typename std::decay_t<T>::value_type& >;
template<class T> //T = std::map
using get_mapped_type = std::conditional_t< std::is_const< std::remove_reference_t<T> >::value,
    const typename std::decay_t<T>::mapped_type&, typename std::decay_t<T>::mapped_type& >;

template<class T>//T = MPIVector
using get_data_type = std::conditional_t< std::is_const< std::remove_reference_t<T> >::value,
    const typename std::decay_t<T>::container_type&, typename std::decay_t<T>::container_type& >;

template<class T>
inline get_element_type<T> do_get_vector_element( T&& v, unsigned i, RecursiveVectorTag)//-> decltype(v[i]){
{
    return v[i];
}
template<class T, class Key>
inline get_mapped_type<T> do_get_vector_element( T&& v, const Key& key, StdMapTag)//-> decltype(v[i]){
{
    return v.at(key);
}
template<class T, class Key>
inline T&& do_get_vector_element( T&& v, const Key&, AnyScalarTag){
    return std::forward<T>(v);
}

template<class T>
inline get_data_type<T> do_get_data( T&& v, MPIVectorTag)//-> decltype(v.data())
{
    return v.data();
}
template<class T>
inline T&& do_get_data( T&& v, AnyScalarTag){
    return std::forward<T>(v);
}

template<class T>
inline get_pointer_type<T> do_get_pointer_or_reference( T&& v, AnyVectorTag)// -> decltype(thrust::raw_pointer_cast(v.data())) //nvcc-7.5 does not like decltype in this connection
{
    return thrust::raw_pointer_cast(v.data());
}
template<class T>
inline T&& do_get_pointer_or_reference( T&& v, AnyScalarTag){
    return std::forward<T>(v);
}

///@endcond

//template<class T>
//inline std::conditional_t<std::is_base_of<AnyScalarTag, get_tensor_category<T>>::value, T&&, get_element_type<T> > get_vector_element( T&& v, unsigned i )// -> decltype( do_get_vector_element( std::forward<T>(v), i, get_tensor_category<T>()) )
//{
//    return do_get_vector_element( std::forward<T>(v), i, get_tensor_category<T>());
//}
//
//template<class T>
//inline std::conditional_t<std::is_base_of<AnyScalarTag, get_tensor_category<T>>::value, T, get_data_type<T> > get_data( T&& v)//-> decltype(do_get_data( std::forward<T>(v), get_tensor_category<T>() ))
//{
//    return do_get_data( std::forward<T>(v), get_tensor_category<T>());
//}
//
//template<class T>
//std::conditional<std::is_base_of<AnyScalarTag, get_tensor_category<T>>::value, T, get_pointer_type<T> > get_pointer_or_reference( T&& v )// -> decltype( do_get_pointer_or_reference( std::forward<T>(v), get_tensor_category<T>()))
//{
//    return do_get_pointer_or_reference( std::forward<T>(v), get_tensor_category<T>());
//}

}//namespace dg

#endif //_DG_VECTOR_TRAITS_
