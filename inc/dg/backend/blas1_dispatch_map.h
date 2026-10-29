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
            std::forward<Params>(ps)...));
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
        dg::assign( el.second, to[el.first], std::forward<Params>(ps)...);
}
} //namespace detail

namespace blas1
{

namespace detail
{

template< class T, size_t N, class Functor, class ContainerType, class ...ContainerTypes>
inline void doDot_fpe( StdMapTag, int* status, std::array<T,N>& fpe, Functor f,
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
        doDot_fpe( status, tmp, f,
            do_get_vector_element(x, el.first, get_tensor_category<ContainerType>()),
            do_get_vector_element(xs, el.first, get_tensor_category<ContainerTypes>())...);
        }catch( std::out_of_range& err)
        {
            throw dg::Error( dg::Message(_ping_)<<"Wrong key '"<<el.first<<"' in blas1::vdot "<<err.what());
        }
        for (unsigned k = 0; k < N; ++k)
            exblas::cpu::Accumulate( tmp[k], fpe, status);
    }
}

template< class Vector1, class Vector2>
inline std::vector<int64_t> doDot_superacc( int* status, const Vector1& x1, const Vector2& x2, StdMapTag)
{
    //find out which one is the RecursiveVector and determine size
    constexpr unsigned vector_idx = find_if_v<dg::is_not_scalar, Vector1, Vector1, Vector2>::value;
    std::vector<int64_t> acc( exblas::BIN_COUNT, (int64_t)0);
    int i=0;
    for( const auto& el : get_idx<vector_idx>(x1,x2))
    {
        try{
        std::vector<int64_t> temp = doDot_superacc(
            status,
            do_get_vector_element(x1,el.first,get_tensor_category<Vector1>()),
            do_get_vector_element(x2,el.first,get_tensor_category<Vector2>()));
        int imin = exblas::IMIN, imax = exblas::IMAX;
        exblas::cpu::Normalize( &(temp[0]), imin, imax);
        for( int k=exblas::IMIN; k<=exblas::IMAX; k++)
            acc[k] += temp[k];
        if( (i+1)%128 == 0)
        {
            imin = exblas::IMIN, imax = exblas::IMAX;
            exblas::cpu::Normalize( &(acc[0]), imin, imax);
        }
        i++;
        }catch( std::out_of_range& err)
        {
            throw dg::Error( dg::Message(_ping_)<<"Wrong key '"<<el.first<<"' in blas1::dot "<<err.what());
        }
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
    for( const auto& el : x)
        result = op( result, doReduce( get_tensor_category<decltype(
            el.second)>(), el.second, zero, op, unary_op));
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
        dg::blas1::subroutine( f, do_get_vector_element(std::forward<container>(x),el.first,get_tensor_category<container>()), do_get_vector_element(std::forward<Containers>(xs),el.first,get_tensor_category<Containers>())...);
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
