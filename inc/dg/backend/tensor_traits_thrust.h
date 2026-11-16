#pragma once

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/complex.h>
#include "vector_categories.h"
#include "tensor_traits.h"
#include "tensor_traits_scalar.h"

namespace dg
{
///@addtogroup traits
///@{

///@brief prototypical Shared Vector with Serial Tag
template<class T>
struct TensorTraits<thrust::host_vector<T>> //, std::enable_if_t< dg::is_scalar<T>::value>>
{
    using value_type        = T; //!< see \c dg::TensorTraits
    using tensor_category   = ThrustVectorTag; //!< see \c dg::TensorTraits
    using execution_policy  = SerialTag; //!< see \c dg::TensorTraits
};
//// MW Thrust should not be Recursive because of user-provided value_types
//// or std::complex or similar
//template<class T>
//struct TensorTraits<thrust::host_vector<T>,
//    std::enable_if_t< !dg::is_scalar<T>::value>>
//{
//    using value_type        = get_value_type<T>;
//    using tensor_category   = RecursiveVectorTag;
//    using execution_policy  = get_execution_policy<T>;
//};

///@brief prototypical Shared Vector. Its execution policy depends on the \c THRUST_DEVICE_SYSTEM macro and can be \c CudaTag, \c OmpTag or \c SerialTag
template<class T>
struct TensorTraits<thrust::device_vector<T> >
{
    using value_type        = T; //!< see \c dg::TensorTraits
    using tensor_category   = ThrustVectorTag; //!< see \c dg::TensorTraits
#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_CUDA
    using execution_policy  = CudaTag ; //!< Depending on the \c THRUST_DEVICE_SYSTEM macro this is either \c CudaTag, \c OmpTag or \c SerialTag. see \c dg::TensorTraits
#elif THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_OMP
    using execution_policy  = OmpTag ; //!< see \c dg::TensorTraits
#else
    using execution_policy  = SerialTag ; //!< see \c dg::TensorTraits
#endif
};

///@brief Enable \c thrust::complex as a Scalar
template<class T>
struct TensorTraits<thrust::complex<T>>
{
    using value_type        = thrust::complex<T>; //!< see \c dg::TensorTraits
    using tensor_category   = ComplexTag; //!< see \c dg::TensorTraits
    using execution_policy  = AnyPolicyTag; //!< see \c dg::TensorTraits
};
///@}
///@cond
//thrust::cpp is an alias for thrust::system::cpp
//tag is a class deriving from thrust::execution_policy<tag>
//raw pointers have the std::random_access_iterator_tag as iterator_category in thrust
template<class Tag>
struct ThrustTag { };
template <>
struct ThrustTag<SerialTag>
{
    using thrust_tag = thrust::cpp::tag;
};
#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_CUDA
template <>
struct ThrustTag<CudaTag>
{
    using thrust_tag = thrust::cuda::tag;
};
#elif THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_OMP
template <>
struct ThrustTag<OmpTag>
{
    using thrust_tag = thrust::omp::tag;
};
#endif
template<class Vector>
using get_thrust_tag = typename ThrustTag<get_execution_policy<Vector>>::thrust_tag;

///@endcond
} //namespace dg
