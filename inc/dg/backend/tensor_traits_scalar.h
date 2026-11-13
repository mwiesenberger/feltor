#pragma once

#include "scalar_categories.h"
#include "tensor_traits.h"

namespace dg
{
///@addtogroup traits
///@{

/// Enable double and float as a floating point
template<class T>
struct TensorTraits<T, std::enable_if_t< std::is_floating_point_v<T>>>
{
    using value_type        = T; //!< see \c dg::TensorTraits
    using tensor_category   = FloatingPointTag; //!< see \c dg::TensorTraits
    using execution_policy  = AnyPolicyTag; //!< see \c dg::TensorTraits
};
/// Enable integers and anything promotable to integer (such as bool and char) as integral
template<class T>
struct TensorTraits<T, std::enable_if_t< std::is_integral_v<T>>>
{
    using value_type        = T; //!< see \c dg::TensorTraits
    using tensor_category   = IntegralTag; //!< see \c dg::TensorTraits
    using execution_policy  = AnyPolicyTag; //!< see \c dg::TensorTraits
};

///@}
} //namespace dg
