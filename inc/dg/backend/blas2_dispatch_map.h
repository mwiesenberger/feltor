#pragma once

#include "blas1_dispatch_map.h"
#include "matrix_categories.h"
#include "../subroutines.h" // dg::Product

///@cond
namespace dg
{
namespace blas2
{
namespace detail
{
template< class ContainerType1, class MatrixType, class ContainerType2>
inline std::vector<int64_t> doDot_superacc( int* status, const ContainerType1& x, const MatrixType& m, const ContainerType2& y);

// x and y are maps and m is not: sum dot( x[key], m, y[key])
template< class Vector1, class Matrix, class Vector2>
inline std::vector<int64_t> doDot_superacc( StdMapTag, int* status,
    const Vector1& x, const Matrix& m, const Vector2& y)
{
    constexpr unsigned vector_idx = find_if_v<dg::is_not_scalar, Vector1, Vector1, Vector2>::value;
    constexpr unsigned K = exblas::ValueTraits<get_result_type<dg::Product, Vector1, Matrix, Vector2>>::K;
    std::vector<int64_t> acc( K*exblas::BIN_COUNT, (int64_t)0);
    *status = 0;
    for( const auto& el : get_idx<vector_idx>(x,y))
    {
        int status_i = 0;
        std::vector<int64_t> temp;
        try{
        temp = doDot_superacc( &status_i,
            do_get_vector_element(x, el.first, get_tensor_category<Vector1>()), m,
            do_get_vector_element(y, el.first, get_tensor_category<Vector2>()));
        }catch( std::out_of_range& err)
        {
            throw dg::Error( dg::Message(_ping_)<<"Wrong key '"<<el.first<<"' in blas2::dot "<<err.what());
        }
        exblas::cpu::AddSuperaccs( acc.data(), temp.data(), K);
        if( status_i != 0)
            *status = status_i;
    }
    return acc;
}
} //namespace detail
} //namespace blas2
} //namespace dg
///@endcond
