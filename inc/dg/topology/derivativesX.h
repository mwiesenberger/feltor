#pragma once

#include "gridX.h"
#include "dxX.h"
#include "../blas.h"

/*! @file
  @brief Convenience functions to create 2D derivatives on X-point topology
  */
namespace dg{

/**
 * @brief A matrix composed of one or two matrices that write to disjoint parts of the output
 *
 * Used for the derivatives on X-point grids (see \c dg::create::derivative
 * and \c dg::create::jump), where the inner and the outer part of the domain
 * need different matrices: applying \c Composite applies \c m1 and then, if
 * \c dual, \c m2 to the same vectors. The two matrices must therefore write
 * to disjoint parts of the output (e.g. via \c dg::EllSparseBlockMat::right_range).
 * @tparam Matrix the type of the two matrices
 * @ingroup creation
 */
template<class Matrix>
struct Composite
{
    ///@brief Empty matrices, not dual
    Composite( ):m1(), m2(), dual(false){ }
    ///@brief Copy from a composite of another matrix type (\c Matrix must be constructible from \c Matrix2)
    ///@param src the composite to copy
    template<class Matrix2>
    Composite( const Composite<Matrix2>& src):m1(src.m1), m2(src.m2), dual(src.dual){}
    ///@brief A single matrix (not dual)
    ///@param m the matrix
    Composite( const Matrix& m):m1(m), m2(m), dual(false){ }
    ///@brief Two matrices (dual)
    ///@param m1 applied first
    ///@param m2 applied second; must write to a part of the output disjoint from \c m1
    Composite( const Matrix& m1, const Matrix& m2):m1(m1), m2(m2), dual(true){ }
    ///@brief Assign from a composite of another matrix type
    ///@param src the composite to copy
    ///@return \c *this
    template<class Matrix2>
    Composite& operator=( const Composite<Matrix2>& src){ Composite c(src);
        *this = c; return *this;}
    ///@brief \f$ v_2 = M v_1\f$ (applies \c m1 and, if \c dual, \c m2)
    ///@param v1 input
    ///@param v2 output
    template< class ContainerType1, class ContainerType2>
    void symv( const ContainerType1& v1, ContainerType2& v2) const
    {
        dg::blas2::symv( m1, v1, v2); //computes first part
        if( dual)
            dg::blas2::symv( m2, v1, v2); //computes second part
    }
    ///@brief \f$ v_2 = \alpha M v_1 + \beta v_2\f$ (applies \c m1 and, if \c dual, \c m2)
    ///@param alpha scalar
    ///@param v1 input
    ///@param beta scalar
    ///@param v2 output
    template< class ContainerType>
    void symv( get_value_type<ContainerType> alpha, const  ContainerType& v1, get_value_type<ContainerType> beta, ContainerType& v2) const
    {
        dg::blas2::symv( alpha, m1, v1, beta, v2); //computes first part
        if( dual)
            dg::blas2::symv( alpha, m2, v1, beta, v2); //computes second part
    }
    ///@brief Display the matrices to a stream
    ///@param os the stream
    void display( std::ostream& os = std::cout) const
    {
        if( dual)
        {
            os << " dual matrix: \n";
            os << " INNER MATRIX\n";
            m1.display( os);
            os << " OUTER MATRIX\n";
            m2.display( os);
        }
        else
        {
            os << "single matrix: \n";
            m1.display(os);
        }
    }
    Matrix m1; //!< the first matrix
    Matrix m2; //!< the second matrix (only used if \c dual)
    bool dual; //!< true if \c m2 is applied
};
///@cond
template <class Matrix>
struct TensorTraits<Composite<Matrix> >
{
    using value_type = get_value_type<Matrix>;
    using tensor_category = SelfMadeMatrixTag;
};
///@endcond


/**
 * @brief Contains functions used for matrix creation
 */
namespace create{

///@addtogroup creation
///@{

//dx, dy, jumpX, jumpY

/**
 * @brief Create 2d derivative in x-direction
 *
 * @param coord Either 0 (for x derivative) or 1 (for y derivative)
 * @param g The grid on which to create dx
 * @param bc The boundary condition
 * @param dir The direction of the first derivative
 *
 * @return A host matrix
 */
template<class real_type>
Composite<EllSparseBlockMat<real_type, thrust::host_vector> > derivative( unsigned coord, const aRealTopologyX2d<real_type>& g, dg::bc bc, direction dir = centered)
{
    if( coord == 0)
    {
        EllSparseBlockMat<real_type, thrust::host_vector>  dx;
        dx = dx_normed( g.n(), g.Nx(), g.hx(), bc, dir);
        dx.set_left_size( g.n()*g.Ny());
        return dx;
    }
    EllSparseBlockMat<real_type, thrust::host_vector>  dy_inner, dy_outer;
    RealGridX1d<real_type> g1d_inner( g.y0(), g.y1(), g.fy(), g.n(), g.Ny(), bc);
    dy_inner = dx( g1d_inner, bc, dir);
    dy_outer = dx_normed( g.n(), g.Ny(), g.hy(), bc, dir );
    dy_inner.right_size = g.n()*g.Nx();
    dy_inner.right_range[0] = 0;
    dy_inner.right_range[1] = g.n()*g.inner_Nx();
    dy_outer.right_range[0] = g.n()*g.inner_Nx();
    dy_outer.right_range[1] = g.n()*g.Nx();
    dy_outer.right_size = g.n()*g.Nx();

    Composite<EllSparseBlockMat<real_type, thrust::host_vector> > c( dy_inner, dy_outer);
    return c;
    // TODO throw on coord > 1 ?
}

/**
 * @brief Matrix that contains 2d jump terms in X direction
 *
 * @param coord Either 0 (for x derivative) or 1 (for y derivative)
 * @param g grid
 * @param bc boundary condition in x
 *
 * @return A host matrix
 */
template<class real_type>
Composite<EllSparseBlockMat<real_type, thrust::host_vector> > jump( unsigned coord, const aRealTopologyX2d<real_type>& g, bc bc)
{
    if( coord == 0)
    {
        EllSparseBlockMat<real_type, thrust::host_vector>  jx;
        jx = jump( g.n(), g.Nx(), g.hx(), bc);
        jx.set_left_size( g.n()*g.Ny());
        return jx;
    }
    EllSparseBlockMat<real_type, thrust::host_vector>  jy_inner, jy_outer;
    RealGridX1d<real_type> g1d_inner( g.y0(), g.y1(), g.fy(), g.n(), g.Ny(), bc);
    jy_inner = jump( g1d_inner, bc);
    jy_outer = jump( g.n(), g.Ny(), g.hy(), bc);
    jy_inner.right_size = g.n()*g.Nx();
    jy_inner.right_range[0] = 0;
    jy_inner.right_range[1] = g.n()*g.inner_Nx();
    jy_outer.right_range[0] = g.n()*g.inner_Nx();
    jy_outer.right_range[1] = g.n()*g.Nx();
    jy_outer.right_size = g.n()*g.Nx();

    Composite<EllSparseBlockMat<real_type, thrust::host_vector> > c( jy_inner, jy_outer);
    return c;
    // TODO throw on coord > 1 ?
}

///////////////////////////////////////////3D VERSIONS//////////////////////
/**
 * @brief Matrix that contains jump terms in X direction in 3D
 *
 * @param coord Either 0 (for x derivative) or 1 (for y derivative) or 2 (for z derivative)
 * @param g The 3D grid
 * @param bc boundary condition
 *
 * @return A host matrix
 */
template<class real_type>
Composite<EllSparseBlockMat<real_type, thrust::host_vector> > jump( unsigned coord, const aRealTopologyX3d<real_type>& g, bc bc)
{
    if( coord == 0)
    {
        EllSparseBlockMat<real_type, thrust::host_vector>  jx;
        jx = jump( g.n(), g.Nx(), g.hx(), bc);
        jx.set_left_size( g.n()*g.Ny()*g.Nz());
        return jx;
    }
    else if ( coord == 1)
    {
        EllSparseBlockMat<real_type, thrust::host_vector>  jy_inner, jy_outer;
        RealGridX1d<real_type> g1d_inner( g.y0(), g.y1(), g.fy(), g.n(), g.Ny(), bc);
        jy_inner = jump( g1d_inner, bc);
        jy_outer = jump( g.n(), g.Ny(), g.hy(), bc);
        jy_inner.right_size = g.n()*g.Nx();
        jy_inner.right_range[0] = 0;
        jy_inner.right_range[1] = g.n()*g.inner_Nx();
        jy_outer.right_range[0] = g.n()*g.inner_Nx();
        jy_outer.right_range[1] = g.n()*g.Nx();
        jy_outer.right_size = g.n()*g.Nx();
        jy_inner.left_size = g.Nz();
        jy_outer.left_size = g.Nz();

        Composite<EllSparseBlockMat<real_type, thrust::host_vector> > c( jy_inner, jy_outer);
        return c;
    }
    EllSparseBlockMat<real_type, thrust::host_vector>  jz;
    jz = jump( 1, g.Nz(), g.hz(), bc);
    jz.set_right_size( g.n()*g.Nx()*g.n()*g.Ny());
    return jz;
    // TODO throw on coord > 2 ?
}


/**
 * @brief Create 3d derivative in x-direction
 *
 * @param coord Either 0 (for x derivative) or 1 (for y derivative) or 2 (for z derivative)
 * @param g The grid on which to create dx
 * @param bc The boundary condition
 * @param dir The direction of the first derivative
 *
 * @return A host matrix
 */
template<class real_type>
Composite<EllSparseBlockMat<real_type, thrust::host_vector> > derivative( unsigned coord, const aRealTopologyX3d<real_type>& g, bc bc, direction dir = centered)
{
    if( coord == 0)
    {
        EllSparseBlockMat<real_type, thrust::host_vector>  dx;
        dx = dx_normed( g.n(), g.Nx(), g.hx(), bc, dir);
        dx.set_left_size( g.n()*g.Ny()*g.Nz());
        return dx;
    }
    else if( coord == 1)
    {
        EllSparseBlockMat<real_type, thrust::host_vector>  dy_inner, dy_outer;
        RealGridX1d<real_type> g1d_inner( g.y0(), g.y1(), g.fy(), g.n(), g.Ny(), bc);
        dy_inner = dx( g1d_inner, bc, dir);
        dy_outer = dx_normed( g.n(), g.Ny(), g.hy(), bc, dir );
        dy_inner.right_size = g.n()*g.Nx();
        dy_inner.right_range[0] = 0;
        dy_inner.right_range[1] = g.n()*g.inner_Nx();
        dy_outer.right_range[0] = g.n()*g.inner_Nx();
        dy_outer.right_range[1] = g.n()*g.Nx();
        dy_outer.right_size = g.n()*g.Nx();
        dy_inner.left_size = g.Nz();
        dy_outer.left_size = g.Nz();

        Composite<EllSparseBlockMat<real_type, thrust::host_vector> > c( dy_inner, dy_outer);
        return c;
    }
    EllSparseBlockMat<real_type, thrust::host_vector>  dz;
    dz = dx_normed( 1, g.Nz(), g.hz(), bc, dir);
    dz.set_right_size( g.n()*g.n()*g.Nx()*g.Ny());
    return dz;
    // TODO throw on coord > 2 ?
}



///@}

} //namespace create

} //namespace dg

