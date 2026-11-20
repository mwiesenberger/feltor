#ifndef _DG_BLAS_PRECONDITIONER_
#define _DG_BLAS_PRECONDITIONER_

#include "tensor_traits.h"
#include "matrix_categories.h"
#include "vector_categories.h"
#include "blas1_dispatch_shared.h"
#include "blas1_dispatch_vector.h"

#include "blas2_stencil.h"

///@cond
namespace dg{

template<class from_ContainerType, class ContainerType, class ...Params>
inline void assign( const from_ContainerType& from, ContainerType& to, Params&& ... ps);

namespace blas1{
template< class ContainerType, class ContainerType1, class ContainerType2, class value_type, class value_type1>
inline void pointwiseDot( value_type alpha, const ContainerType1& x1, const ContainerType2& x2, value_type1 beta, ContainerType& y);
template< class ContainerType, class ContainerType1, class ContainerType2>
inline void pointwiseDot( const ContainerType1& x1, const ContainerType2& x2, ContainerType& y);
}//namespace blas1

namespace blas2{
//forward declare blas2 functions
template< class MatrixType, class ContainerType1, class ContainerType2>
inline void symv( get_value_type<ContainerType1> alpha,
                  MatrixType&& M,
                  const ContainerType1& x,
                  get_value_type<ContainerType1> beta,
                  ContainerType2& y);
template< class MatrixType, class ContainerType1, class ContainerType2>
inline void symv( MatrixType&& M,
                  const ContainerType1& x,
                  ContainerType2& y);
namespace detail{

//thrust vector preconditioner
template< class Vector1, class Vector2>
void doTransfer( const Vector1& in, Vector2& out, AnyVectorTag, AnyVectorTag)
{
    dg::assign(in,out);
}

template< class Stencil, class ContainerType, class ...ContainerTypes>
inline void doParallelFor( SharedVectorTag, Stencil f, unsigned N, ContainerType&& x, ContainerTypes&&... xs)
{
    // a copy of doSubroutine ...

    using vector_type = find_if_t<dg::is_not_scalar_has_not_any_policy, get_value_type<ContainerType>, ContainerType, ContainerTypes...>;
    using execution_policy = get_execution_policy<vector_type>;
    static_assert((
            dg::has_any_or_same_policy<ContainerType, execution_policy>::value &&
      ... &&dg::has_any_or_same_policy<ContainerTypes, execution_policy>::value),
        "All ContainerType types must have compatible execution policies (AnyPolicy or Same)!");
    doParallelFor_dispatch(
            get_execution_policy<vector_type>(),
            N,
            f,
            do_get_pointer_or_reference(std::forward<ContainerType>(x),get_tensor_category<ContainerType>()) ,
            do_get_pointer_or_reference(std::forward<ContainerTypes>(xs),get_tensor_category<ContainerTypes>()) ...
            );
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              get_value_type<Vector1> alpha,
              Matrix&& m,
              const Vector1& x,
              get_value_type<Vector1> beta,
              Vector2& y,
              SharedVectorTag, SharedVectorTag)
{
    dg::blas1::pointwiseDot( alpha, std::forward<Matrix>(m), x, beta, y);
}
template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              Matrix&& m,
              const Vector1& x,
              Vector2& y,
              SharedVectorTag, SharedVectorTag)
{
    dg::blas1::pointwiseDot( std::forward<Matrix>(m), x, y);
}
template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              get_value_type<Vector1> alpha,
              Matrix&& m,
              const Vector1& x,
              get_value_type<Vector1> beta,
              Vector2& y,
              SharedVectorTag, RecursiveVectorTag)
{
    for(unsigned i=0; i<x.size(); i++)
        dg::blas2::symv( alpha, m, do_get_vector_element(x,i,get_tensor_category<Vector1>()), beta, do_get_vector_element(y,i,get_tensor_category<Vector2>()));
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              Matrix&& m,
              const Vector1& x,
              Vector2& y,
              SharedVectorTag, RecursiveVectorTag)
{
    for(unsigned i=0; i<y.size(); i++)
        dg::blas2::symv( m, do_get_vector_element(x,i,get_tensor_category<Vector1>()), do_get_vector_element(y,i,get_tensor_category<Vector2>()));
}
template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              get_value_type<Vector1> alpha,
              Matrix&& m,
              const Vector1& x,
              get_value_type<Vector1> beta,
              Vector2& y,
              SharedVectorTag, StdMapTag)
{
    for( const auto& el : y)
        dg::blas2::symv( alpha, m, do_get_vector_element(x,el.first,get_tensor_category<Vector1>()), beta, do_get_vector_element(y,el.first,get_tensor_category<Vector2>()));
}

template< class Matrix, class Vector1, class Vector2>
inline void doSymv(
              Matrix&& m,
              const Vector1& x,
              Vector2& y,
              SharedVectorTag, StdMapTag)
{
    for( const auto& el : y)
        dg::blas2::symv( m, do_get_vector_element(x,el.first,get_tensor_category<Vector1>()), do_get_vector_element(y,el.first,get_tensor_category<Vector2>()));
}


}//namespace detail
} //namespace blas2
} //namespace dg
///@endcond
#endif //_DG_BLAS_PRECONDITIONER_
