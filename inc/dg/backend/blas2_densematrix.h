#pragma once

#include <tuple>
#include <utility>
#include "densematrix.h"

namespace dg{
namespace blas2{
namespace detail{

// one blas1::evaluate over the columns offset, ..., offset+sizeof...(I)-1
template< class Matrix, class Vector1, class Vector2, std::size_t ...I>
inline void dense_pass( get_value_type<Vector1> alpha, const Matrix& m,
    const Vector1& x, unsigned offset, get_value_type<Vector1> beta, Vector2& y,
    std::index_sequence<I...>)
{
    using value_type = get_value_type<Vector1>;
    using container = typename std::decay_t<Matrix>::container_type;
    std::apply( [&]( const auto& ... args){
        dg::blas1::evaluate( y, dg::Axpby(alpha, beta), dg::PairSum(), args...);
    }, std::tuple_cat( std::tuple<value_type, const container&>(
        x[offset+I], *m.get()[offset+I])...));
}

// s.a. the detail::gemm algorithm in runge_kutta.h
template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              get_value_type<Vector1> alpha,
              Matrix&& m,
              const Vector1& x,
              get_value_type<Vector1> beta,
              Vector2& y,
              DenseMatrixTag)
{
    static_assert( std::is_same<get_execution_policy<Matrix>,
                                get_execution_policy<Vector2>>::value,
                "Dense Matrix and output Vector type must have same execution policy");
    static_assert( std::is_same<get_tensor_category<
                            typename std::decay_t<Matrix>::container_type>,
                        get_tensor_category<Vector2>>::value,
                "Dense Matrix and output Vector type must have same data layout");
    const unsigned size = x.size();
    if( size != m.num_cols()) {
        throw Error( Message(_ping_)<<"x has the wrong size "<<x.size()<<" Number of columns is "<<m.num_cols());
    }
    if( size == 0)
    {
        dg::blas1::scal( y, beta);
        return;
    }
    // Sum 8 columns per pass, then the remaining columns in one pass
    // (every pass reads and writes y)
    unsigned i=0;
    for( ; i+8 <= size; i+=8)
        dense_pass( alpha, m, x, i, i == 0 ? beta : 1., y, std::make_index_sequence<8>());
    const get_value_type<Vector1> b = i == 0 ? beta : 1.;
    switch( size - i)
    {
        case 1: dg::blas1::axpby( alpha*x[i], *m.get()[i], b, y); break;
        case 2: dense_pass( alpha, m, x, i, b, y, std::make_index_sequence<2>()); break;
        case 3: dense_pass( alpha, m, x, i, b, y, std::make_index_sequence<3>()); break;
        case 4: dense_pass( alpha, m, x, i, b, y, std::make_index_sequence<4>()); break;
        case 5: dense_pass( alpha, m, x, i, b, y, std::make_index_sequence<5>()); break;
        case 6: dense_pass( alpha, m, x, i, b, y, std::make_index_sequence<6>()); break;
        case 7: dense_pass( alpha, m, x, i, b, y, std::make_index_sequence<7>()); break;
        default: break;
    }
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              Matrix&& m,
              const Vector1& x,
              Vector2& y,
              DenseMatrixTag)
{
    doSymv( 1., std::forward<Matrix>(m), x, 0., y, DenseMatrixTag());
}

} //namespace detail
} //namespace blas2
} //namespace dg
