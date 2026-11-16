#pragma once

#include <array>
#include <vector>
#include <map>
#include <complex>

#include "vector_categories.h"
#include "tensor_traits.h"
#include "tensor_traits_scalar.h"
#include "predicate.h"

namespace dg
{
///@addtogroup traits
///@{

///@brief Prototypical Recursive Vector (unless vector of scalars)
template<class T>
struct TensorTraits<std::vector<T>,
    std::enable_if_t< !dg::is_scalar<T>::value>>
{
    using value_type        = get_value_type<T>; //!< see \c dg::TensorTraits
    using tensor_category   = RecursiveVectorTag; //!< see \c dg::TensorTraits
    using execution_policy  = get_execution_policy<T>; //!< see \c dg::TensorTraits
};
///@brief Prototypical Shared Vector (vector of scalars) with Serial Tag
template<class T>
struct TensorTraits<std::vector<T>,
    std::enable_if_t< dg::is_scalar<T>::value>>
{
    using value_type        = T; //!< see \c dg::TensorTraits
    using tensor_category   = ThrustVectorTag; //!< see \c dg::TensorTraits
    using execution_policy  = SerialTag; //!< see \c dg::TensorTraits
};

///@brief Prototypical Recursive Vector of fixed size (unless array of scalars)
template<class T, std::size_t N>
struct TensorTraits<std::array<T, N>,
    std::enable_if_t< !dg::is_scalar<T>::value>>
{
    using value_type        = get_value_type<T>; //!< see \c dg::TensorTraits
    using tensor_category   = ArrayVectorTag; //!< see \c dg::TensorTraits
    using execution_policy  = get_execution_policy<T>; //!< see \c dg::TensorTraits
};

///@brief Array of scalars: a Shared Vector of fixed size with Serial Tag
template<class T, std::size_t N>
struct TensorTraits<std::array<T, N>,
    std::enable_if_t< dg::is_scalar<T>::value>>
{
    using value_type        = T; //!< see \c dg::TensorTraits
    using tensor_category   = ArrayScalarTag; //!< see \c dg::TensorTraits
    using execution_policy  = SerialTag; //!< see \c dg::TensorTraits
};

/// Behaves like a RecursiveVector
template<class Key, class T>
struct TensorTraits<std::map<Key,T>>
{
    using value_type        = get_value_type<T>; //!< see \c dg::TensorTraits
    using tensor_category   = StdMapTag; //!< see \c dg::TensorTraits
    using execution_policy  = get_execution_policy<T>; //!< see \c dg::TensorTraits
};

/// Enable std::complex as a Scalar
template<class T>
struct TensorTraits<std::complex<T>>
{
    using value_type        = std::complex<T>; //!< see \c dg::TensorTraits
    using tensor_category   = ComplexTag; //!< see \c dg::TensorTraits
    using execution_policy  = AnyPolicyTag; //!< see \c dg::TensorTraits
};
///@}
} //namespace dg
