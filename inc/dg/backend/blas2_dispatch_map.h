#pragma once

#include "blas1_dispatch_map.h"
#include "matrix_categories.h"

///@cond
namespace dg
{
namespace blas2
{
namespace detail
{

template< class ContainerType1, class MatrixType, class ContainerType2>
inline std::vector<int64_t> doDot_superacc( int* status, const ContainerType1& x,
    const MatrixType& m, const ContainerType2& y);

// The element of m that belongs to key: if m is a map its element, else m itself
template<class Matrix, class Key>
const auto& get_map_matrix_element( const Matrix& m, const Key& key)
{
    if constexpr( std::is_base_of_v<StdMapTag, get_tensor_category<Matrix>>)
        return m.at( key);
    else
        return m;
}

// MatrixCategory is deduced from the tensor_category of the matrix (any).
// A fixed base like AnyMatrixTag would be ambiguous with e.g. the
// ( SharedVectorTag, RecursiveVectorTag) overload; deduced, it matches exactly
template< class Vector1, class Matrix, class Vector2, class MatrixCategory>
inline std::vector<int64_t> doDot_superacc( int* status, const Vector1& x,
    const Matrix& m, const Vector2& y, MatrixCategory, StdMapTag)
{
    //find out which one is the map and determine the keys
    constexpr unsigned vector_idx = find_if_v<dg::is_not_scalar, Vector1, Vector1, Vector2>::value;
    std::vector<int64_t> acc( exblas::BIN_COUNT, (int64_t)0);
    int i=0;
    for( const auto& el : get_idx<vector_idx>(x,y))
    {
        try{
        std::vector<int64_t> temp = doDot_superacc( status,
            do_get_vector_element(x,el.first,get_tensor_category<Vector1>()),
            get_map_matrix_element( m, el.first),
            do_get_vector_element(y,el.first,get_tensor_category<Vector2>()));
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
            throw dg::Error( dg::Message(_ping_)<<"Wrong key '"<<el.first<<"' in blas2::dot "<<err.what());
        }
    }
    return acc;
}

} //namespace detail
} //namespace blas2
} //namespace dg
///@endcond
