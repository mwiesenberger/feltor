#pragma once
#include <vector>
#include "exblas/mpi_accumulate.h"
#include "mpi_matrix.h"
#include "blas1_dispatch_mpi.h"
#include "blas2_dispatch_shared.h"

///@cond
namespace dg
{
namespace blas2
{
template< class Stencil, class ContainerType, class ...ContainerTypes>
inline void parallel_for( Stencil f, unsigned N, ContainerType&& x, ContainerTypes&&... xs);
namespace detail
{

template< class Matrix1, class Matrix2>
inline void doTransfer( const Matrix1& m1, Matrix2& m2, AnyMatrixTag, MPIMatrixTag)
{
    Matrix2 m(m1);
    m2 = m;
}

//Matrix = MPI Vector
template< class Matrix, class Vector1, class Vector2>
inline void doSymv( Matrix&& m, const Vector1& x, Vector2& y, MPIVectorTag, MPIVectorTag )
{
    dg::blas2::symv( m.data(), do_get_data(x, get_tensor_category<Vector1>()), do_get_data(y, get_tensor_category<Vector2>()));
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv( get_value_type<Vector1> alpha,
                Matrix&& m,
                const Vector1& x,
                get_value_type<Vector1> beta,
                Vector2& y,
                MPIVectorTag,
                MPIVectorTag
                )
{
    dg::blas2::symv( alpha, m.data(), do_get_data(x,get_tensor_category<Vector1>()), beta, do_get_data(y, get_tensor_category<Vector2>()));
}
template< class Matrix, class Vector1, class Vector2>
inline void doSymv( Matrix&& m, const Vector1& x, Vector2& y, MPIVectorTag, RecursiveVectorTag )
{
    for( unsigned i=0; i<y.size(); i++)
        dg::blas2::symv( m, do_get_vector_element(x,i,get_tensor_category<Vector1>()), do_get_vector_element(y,i,get_tensor_category<Vector2>()));
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv( get_value_type<Vector1> alpha,
                Matrix&& m,
                const Vector1& x,
                get_value_type<Vector1> beta,
                Vector2& y,
                MPIVectorTag,
                RecursiveVectorTag
                )
{
    for( unsigned i=0; i<y.size(); i++)
        dg::blas2::symv( alpha, m, do_get_vector_element(x,i,get_tensor_category<Vector1>()), beta, do_get_vector_element(y,i,get_tensor_category<Vector2>()));
}
//Matrix is an MPI matrix
template< class Matrix, class Vector1, class Vector2>
inline void doSymv( Matrix&& m, const Vector1& x, Vector2& y, MPIMatrixTag, MPIVectorTag )
{
    m.symv( x, y);
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv( get_value_type<Vector1> alpha,
                Matrix&& m,
                const Vector1& x,
                get_value_type<Vector1> beta,
                Vector2& y,
                MPIMatrixTag,
                MPIVectorTag
                )
{
    m.symv( alpha, x, beta, y);
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv( Matrix&& m, const Vector1& x, Vector2& y, MPIMatrixTag, RecursiveVectorTag )
{
    for( unsigned i=0; i<y.size(); i++)
        dg::blas2::symv( m, do_get_vector_element(x,i,get_tensor_category<Vector1>()), do_get_vector_element(y,i,get_tensor_category<Vector2>()));
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv( get_value_type<Vector1> alpha,
                Matrix&& m,
                const Vector1& x,
                get_value_type<Vector1> beta,
                Vector2& y,
                MPIMatrixTag,
                RecursiveVectorTag
                )
{
    for( unsigned i=0; i<y.size(); i++)
        dg::blas2::symv( alpha, m, x[i], beta, y[i]);
}
template< class Matrix, class Vector1, class Vector2>
inline void doSymv( Matrix&& m, const Vector1& x, Vector2& y, MPIMatrixTag, StdMapTag )
{
    for( const auto& el : y)
        dg::blas2::symv( m, do_get_vector_element(x,el.first,get_tensor_category<Vector1>()), do_get_vector_element(y,el.first,get_tensor_category<Vector2>()));
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv( get_value_type<Vector1> alpha,
                Matrix&& m,
                const Vector1& x,
                get_value_type<Vector1> beta,
                Vector2& y,
                MPIMatrixTag,
                StdMapTag
                )
{
    for( const auto& el : y)
        dg::blas2::symv( alpha, m, do_get_vector_element(x,el.first,get_tensor_category<Vector1>()), beta, do_get_vector_element(y,el.first,get_tensor_category<Vector2>()));
}

template<class Functor, class Matrix, class Vector1, class Vector2>
inline void doStencil(
                    Functor f,
                    Matrix&& m,
                    const Vector1&x,
                    Vector2& y,
                    MPIMatrixTag,
                    MPIVectorTag  )
{
    m.stencil( f, x,y);
}


} //namespace detail
} //namespace blas2
} //namespace dg
///@endcond
