#pragma once

#include "blas1_dispatch_vector.h"
#include "matrix_categories.h"
#include "../subroutines.h" // dg::Product

///@cond
namespace dg
{
namespace blas2
{
//forward declare blas2 symv functions
template< class MatrixType, class ContainerType1, class ContainerType2>
void symv( MatrixType&& M,
                  const ContainerType1& x,
                  ContainerType2& y);
template< class MatrixType, class ContainerType1, class ContainerType2>
void symv( get_value_type<ContainerType1> alpha,
                  MatrixType&& M,
                  const ContainerType1& x,
                  get_value_type<ContainerType1> beta,
                  ContainerType2& y);
namespace detail
{
template< class ContainerType1, class MatrixType, class ContainerType2>
inline std::vector<int64_t> doDot_superacc( int* status, const ContainerType1& x, const MatrixType& m, const ContainerType2& y);

// x and y are vectors of containers and m is not: sum dot( x[i], m, y[i])
template< class Vector1, class Matrix, class Vector2>
inline std::vector<int64_t> doDot_superacc( RecursiveVectorTag, int* status,
    const Vector1& x, const Matrix& m, const Vector2& y)
{
    constexpr unsigned vector_idx = find_if_v<dg::is_not_scalar, Vector1, Vector1, Vector2>::value;
    auto size = get_idx<vector_idx>(x,y).size();
    constexpr unsigned K = exblas::ValueTraits<get_result_type<dg::Product, Vector1, Matrix, Vector2>>::K;
    std::vector<int64_t> acc( K*exblas::BIN_COUNT, (int64_t)0);
    *status = 0;
    for( unsigned i=0; i<size; i++)
    {
        int status_i = 0;
        std::vector<int64_t> temp = doDot_superacc( &status_i,
            do_get_vector_element(x, i, get_tensor_category<Vector1>()), m,
            do_get_vector_element(y, i, get_tensor_category<Vector2>()));
        exblas::cpu::AddSuperaccs( acc.data(), temp.data(), K);
        if( status_i != 0)
            *status = status_i;
    }
    return acc;
}


//In case the matrix is a RecursiveVector just do a recursive call
template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              get_value_type<Vector1> alpha,
              Matrix&& m,
              const Vector1& x,
              get_value_type<Vector1> beta,
              Vector2& y,
              RecursiveVectorTag)
{
    for( unsigned i=0; i<m.size(); i++)
        dg::blas2::symv( alpha, m[i], do_get_vector_element(x,i,get_tensor_category<Vector1>()), beta, do_get_vector_element(y,i,get_tensor_category<Vector2>()));
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              Matrix&& m,
              const Vector1& x,
              Vector2& y,
              RecursiveVectorTag)
{
    for( unsigned i=0; i<m.size(); i++)
        dg::blas2::symv( m[i], do_get_vector_element(x,i,get_tensor_category<Vector1>()), do_get_vector_element(y,i,get_tensor_category<Vector2>()));
}



} //namespace detail
} //namespace blas1
} //namespace dg
///@endcond
