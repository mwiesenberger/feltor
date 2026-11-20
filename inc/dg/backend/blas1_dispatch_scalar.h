#ifndef _DG_BLAS_SCALAR_
#define _DG_BLAS_SCALAR_

#include "scalar_categories.h"
#include "tensor_traits.h"
#include "predicate.h"

#include "blas1_serial.h"

///@cond
namespace dg
{
namespace detail
{
template< class To, class From, class ...Params>
To doConstruct( const From& from, AnyScalarTag, AnyScalarTag, Params&& ...)
{
    return To( from);
}
template< class From, class To, class ...Params>
void doAssign( const From& from, To& to, AnyScalarTag, AnyScalarTag, Params&& ...)
{
    to = from;
}
}//namespace detail
namespace blas1
{
namespace detail
{
template< class Functor, class ...ContainerTypes>
inline auto doDot_integral( AnyScalarTag, Functor f, const ContainerTypes& ...xs)
{
    return f( xs...);
}

template< class Functor, class ...ContainerTypes>
inline std::vector<int64_t> doDot_superacc( AnyScalarTag, int* status, Functor f,
    const ContainerTypes& ...xs)
{
    //all are scalars: the sum has a single term, so a small FPE suffices
    constexpr unsigned K = exblas::ValueTraits<get_result_type<Functor, ContainerTypes...>>::K;
    std::vector<int64_t> h_superacc(K*exblas::BIN_COUNT);
    exblas::exsum_cpu<2>( 1, &h_superacc[0], status, f, xs...);
    return h_superacc;
}

template< class Subroutine, class ContainerType, class ...ContainerTypes>
inline void doSubroutine( AnyScalarTag, Subroutine f, ContainerType&& x, ContainerTypes&&... xs)
{
    f(x,xs...);
}

template<class T, class ContainerType, class BinaryOp, class UnaryOp>
inline T doReduce( AnyScalarTag, ContainerType x, T init, BinaryOp op, UnaryOp
        unary_op)
{
    init = op( init, unary_op(x));
    return init;
}
template< class BinarySubroutine, class Functor, class ContainerType, class ...ContainerTypes>
inline void doKronecker( AnyScalarTag, ContainerType& y, BinarySubroutine f, Functor g, const ContainerTypes&... xs)
{
    f( g(xs...),y);
}


} //namespace detail
} //namespace blas1
namespace detail
{
template<class ContainerType, class Functor, class ...ContainerTypes>
auto doKronecker( AnyScalarTag, Functor f, const ContainerType& x0, const ContainerTypes& ... xs)
{
    return f(x0,xs...);
}
} //namespace detail
} //namespace dg
///@endcond

#endif //_DG_BLAS_SCALAR_
