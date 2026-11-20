#pragma once
#include <vector>
#include <array>
#include <map>
#include "blas1_dispatch_shared.h"
#include "vector_categories.h"
#include "tensor_traits.h"
#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP

///@cond
namespace dg
{
template<class to_ContainerType, class from_ContainerType, class ...Params>
inline to_ContainerType construct( const from_ContainerType& src, Params&& ...ps);
template<class from_ContainerType, class to_ContainerType, class ...Params>
inline void assign( const from_ContainerType&, to_ContainerType&, Params&& ...ps);

namespace detail
{
template<class To, class From, class ...Params>
To doConstruct( const From& src, StdMapTag, StdMapTag, Params&&...ps )
{
    To t;
    using inner_vector = typename To::mapped_type;
    for( const auto& el : src)
        t.emplace( el.first, dg::construct<inner_vector>(el.second,
            ps...));
    return t;
}
template<class From, class To, class ...Params>
void doAssign( const From& src, To& to, StdMapTag, StdMapTag, Params&&...ps )
{
    // to gets the same keys as src
    for( auto it = to.begin(); it != to.end(); )
    {
        if( src.count( it->first) == 0)
            it = to.erase( it);
        else
            ++it;
    }
    for( const auto& el : src)
        dg::assign( el.second, to[el.first], ps...);
}
} //namespace detail

namespace blas1
{

namespace detail
{

template< class T, size_t N, class Functor, class ContainerType, class ...ContainerTypes>
inline void doDot_fpe( StdMapTag, std::array<T,N>& fpe, Functor f,
    const ContainerType& x, const ContainerTypes& ...xs)
{
    //find out which one is the map and determine the keys
    constexpr unsigned vector_idx = find_if_v<dg::is_not_scalar, ContainerType, ContainerType, ContainerTypes...>::value;
    //reduce received data (serial execution)
    for( unsigned k=0; k<N; k++)
        fpe[k] = T(0);
    for( const auto& el : get_idx<vector_idx>(x,xs...))
    {
        std::array<T,N> tmp;
        try{
        doDot_fpe( tmp, f,
            do_get_vector_element(x, el.first, get_tensor_category<ContainerType>()),
            do_get_vector_element(xs, el.first, get_tensor_category<ContainerTypes>())...);
        }catch( std::out_of_range& err)
        {
            throw dg::Error( dg::Message(_ping_)<<"Wrong key '"<<el.first<<"' in blas1::vdot "<<err.what());
        }
        for (unsigned k = 0; k < N; ++k)
            exblas::cpu::Accumulate( tmp[k], fpe);
    }
}

template< class Functor, class ...ContainerTypes>
inline std::vector<int64_t> doDot_superacc( StdMapTag, int* status, Functor f,
    const ContainerTypes& ...xs)
{
    //find out which one is the map and determine the keys
    constexpr unsigned vector_idx = find_if_v<dg::is_not_scalar, std::tuple_element_t<0, std::tuple<ContainerTypes...>>, ContainerTypes...>::value;
    constexpr unsigned K = exblas::ValueTraits<get_result_type<Functor, ContainerTypes...>>::K;
    std::vector<int64_t> acc( K*exblas::BIN_COUNT, (int64_t)0);
    *status = 0;
    for( const auto& el : get_idx<vector_idx>(xs...))
    {
        int status_u = 0;
        std::vector<int64_t> temp;
        try{
        temp = doDot_superacc( &status_u, f,
            do_get_vector_element(xs, el.first, get_tensor_category<ContainerTypes>())...);
        }catch( std::out_of_range& err)
        {
            throw dg::Error( dg::Message(_ping_)<<"Wrong key '"<<el.first<<"' in blas1::vdot "<<err.what());
        }
        exblas::cpu::AddSuperaccs( acc.data(), temp.data(), K);
        if( status_u != 0)
            *status = status_u;
    }
    return acc;
}
/////////////////////////////////////////////////////////////////////////////////////

template<class T, class ContainerType, class BinaryOp, class UnaryOp>
inline T doReduce( StdMapTag, const ContainerType& x, T zero, BinaryOp
        op, UnaryOp unary_op)
{
    //reduce sequentially recursively
    // zero must stay the neutral element in the inner reductions (e.g. MPI)
    T result = zero;
    // mapped_type instead of decltype(el.second): nvcc 12.4 fails to compile
    // decltype of the range-for variable in host code
    using inner_vector = typename ContainerType::mapped_type;
    for( const auto& el : x)
        result = op( result, doReduce( get_tensor_category<inner_vector>(),
            el.second, zero, op, unary_op));
    return result;
}

//dispatch
template< class Subroutine, class container, class ...Containers>
inline void doSubroutine( StdMapTag, Subroutine f, container&& x, Containers&&... xs)
{
    constexpr unsigned vector_idx = find_if_v<dg::is_not_scalar, get_value_type<container>, container, Containers...>::value;
    for( const auto& el : get_idx<vector_idx>( std::forward<container>(x), std::forward<Containers>(xs)...))
    {
        try{
        dg::blas1::subroutine( f, do_get_vector_element(x,el.first,get_tensor_category<container>()), do_get_vector_element(xs,el.first,get_tensor_category<Containers>())...);
        }catch( std::out_of_range& err)
        {
            throw dg::Error( dg::Message(_ping_)<<"Wrong key '"<<el.first<<"' in subroutine");
        }
    }
}

} //namespace detail
} //namespace blas1
} //namespace dg
///@endcond
