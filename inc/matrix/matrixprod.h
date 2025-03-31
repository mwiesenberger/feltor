
#pragma once

#include "dg/algorithm.h"
#include "lanczos.h"

namespace dg{
namespace mat{

/**
 * @brief Computation of \f$ \vec x = f(A,\vec d)\vec b\f$ and \f$ \vec x = f(\vec d, A)\vec b\f$
 * where \f$ A \f$ is
 * a positive definite matrix self-adjoint in the weights \f$ W\f$ .
 *
 * The first identity is computed via \f$ \vec x = f(\vec d, A) \vec b = (E_{A} \odot F ) E^T_{A}M^T b\f$
 * where \f$ E_A := V_A E_T \f$ and \f$ F_{ai} := f( d_a, \lambda_i)\f$
 * and \f$ T\f$ and \f$ V_A\f$  are the tridiagonal matrix and vectors that
 * come out of a Lanczos iteration on \f$ A\f$, \f$ W\f$, \f$ \vec b\f$; \f$ \vec d\f$ is a vector.
 *
 * The second identity is computed via \f$ \vec x = f(A, \vec d) \vec b = E_{A} (F^T \odot   E^T_{A}M^T) b\f$
 * where \f$ E_A := V_A E_T \f$ and \f$ F_{ai} := f( d_a, \lambda_i)\f$
 * and \f$ T\f$ and \f$ V_A\f$  are the tridiagonal matrix and vectors that
 * come out of a Lanczos iteration on \f$ A\f$, \f$ W\f$, \f$ \vec b\f$; \f$ \vec d\f$ is a vector
 *
 * @ingroup matrixfunctionapproximation
 * @attention Just as in the Lanczos or PCG methods the matrix \f$ A\f$ needs to be positive-definite (i.e. it won't work for negative definite)
 * @note The \c apply and \c apply_adjoint methods are just abbreviations. If one wants full control, e.g. to reuse a tridiagonalisation one has to manually code:
 *
 * @code{.cpp}
 * double max = dg::blas1::reduce( diag, -1e308, thrust::maximum<double>());
    auto func = dg::mat::make_FuncEigen_Te1( [&](value_type x) {return op( max, x);});
    dg::mat::ProductMatrixFunction<ContainerType> prod( x, 100);
    auto T = prod.lanczos().tridiag( func, A,
                b, weights, eps, nrmb_correction,
                "universal", 1.0, 1);
    prod.compute_vlcl( op, diag, A, T, x, b, prod.lanczos().get_bnorm());
    // or
    prod.compute_vlcl_adjoint( op, A, diag, T, x, b,
                weights, prod.lanczos().get_bnorm());
 * @endcode
 */
template<class ContainerType>
struct ProductMatrixFunction
{
    using container_type = ContainerType;
    using value_type = dg::get_value_type<ContainerType>;
    /// Construct empty
    ProductMatrixFunction() = default;

    /**
     * @brief Allocate memory for the method
     *
     * @param copyable A ContainerType must be copy-constructible from this
     * @param max_iterations Maximum number of iterations to be used
     */
    ProductMatrixFunction( const ContainerType& copyable, unsigned max_iterations)
    {
        m_lanczos.construct( copyable, max_iterations);
        m_v = m_vp = m_vm = m_f = copyable;
    }
    ///@copydoc hide_construct
    template<class ...Params>
    void construct( Params&& ...ps)
    {
        //construct and swap
        *this = ProductMatrixFunction( std::forward<Params>( ps)...);
    }

    ///@copydoc MatrixFunction::set_benchmark(bool,std::string)
    void set_benchmark( bool benchmark, std::string message = "ProductFunction"){
        m_benchmark = benchmark;
        m_message = message;
    }

    /**
     * @brief Compute \f$ \vec x = f(\vec d, A) \vec b = (E_{A} \odot F ) E^T_{A}M^T b\f$
     *
     * This function is equivalent to:
     * @code{.cpp}
        auto func = dg::mat::make_FuncEigen_Te1( [&](value_type x) {return op(1., x);});
        auto T = m_lanczos.tridiag( func, std::forward<MatrixType>(A),
                b, weights, eps, nrmb_correction,
                "universal", 1.0, 2);
        compute_vlcl( op, diag, std::forward<MatrixType>(A), T, x, b,
                    m_lanczos.get_bnorm());
        return T.num_rows;
     * @endcode
     * @note The stopping criterion used on the Lanczos iteration is the
     * universal one applied to \f$ f(1, x) \f$
     * @param x output-vector, contains result on output, ignored on input
     * @param op a  binary Operator representing the product matrix function
     * @param diag the diagonal vector
     * @param A A self-adjoint, positive definit matrix
     * @param b The initial vector that starts orthogonalization
     * @param weights Weights that define the scalar product in which \c A is
     *  self-adjoint and in which the error norm is computed.
     * @param eps relative accuracy of residual in Lanczos iteration
     * @param nrmb_correction the absolute error \c C in units of \c eps to be
     * respected
     * @return The number of Lanczos iterations used
     */
    template<class ContainerType0, class BinaryOp, class ContainerType1,
        class MatrixType, class ContainerType2, class ContainerType3>
    unsigned apply(
            ContainerType0& x,
            BinaryOp op,
            const ContainerType1& diag,
            MatrixType&& A,
            const ContainerType2& b,
            const ContainerType3& weights,
            value_type eps,
            value_type nrmb_correction = 1.)
    {
#ifdef MPI_VERSION
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif //MPI
        dg::Timer t;
        t.tic();
        auto func = make_FuncEigen_Te1( [&](value_type x) {return op(1., x);});
        auto T = m_lanczos.tridiag( func, std::forward<MatrixType>(A),
                b, weights, eps, nrmb_correction,
                "universal", 1.0, 2);
        compute_vlcl( op, diag, std::forward<MatrixType>(A), T, x, b,
                    m_lanczos.get_bnorm());
        t.toc();
        if( m_benchmark)
            DG_RANK0 std::cout << "# `"<<m_message<<"` solve with {"<<T.num_rows<<"} iterations took "<<t.diff()<<"s\n";
        return T.num_rows;
    }
    /**
     * @brief Compute \f$ \vec x = f(A, \vec d) \vec b = E_{A} (F^T \odot   E^T_{A}M^T) b\f$
     *
     * This function is equivalent to:
     * @code{.cpp}
        auto func = make_FuncEigen_Te1( [&](value_type x) {return op( x, 1.);});
        auto T = m_lanczos.tridiag( func, std::forward<MatrixType>(A),
                b, weights, eps, nrmb_correction,
                "universal", 1.0, 2);
        compute_vlcl_adjoint( op, std::forward<MatrixType>(A), diag, T, x, b,
                weights, m_lanczos.get_bnorm());
        return T.num_rows;
     * @endcode
     * @note \f$ f(A, \vec d)\f$ is the adjoint operation to \f$ f( \vec d, A)\f$
     *  since both \f$ \vec d\f$ and \f$ A\f$ are self-adjoint.
     * @note The stopping criterion used on the Lanczos iteration is the
     * universal one applied to \f$ f(x, 1) \f$
     * @param x output-vector, contains result on output, ignored on input
     * @param op a  binary Operator representing the product matrix function
     * @param diag the diagonal vector
     * @param A A self-adjoint, positive definit matrix
     * @attention The order of \c A and \c diag is reversed compared to the
     * \c apply method
     * @param b The initial vector that starts orthogonalization
     * @param weights Weights that define the scalar product in which \c A is
     *  self-adjoint and in which the error norm is computed.
     * @param eps relative accuracy of residual in Lanczos iteration
     * @param nrmb_correction the absolute error \c C in units of \c eps to be
     * respected
     * @return The number of Lanczos iterations used
     */
    template<class ContainerType0, class BinaryOp, class MatrixType,
        class ContainerType1, class ContainerType2, class ContainerType3>
    unsigned apply_adjoint(
            ContainerType0& x,
            BinaryOp op,
            MatrixType&& A,
            const ContainerType1& diag,
            const ContainerType2& b,
            const ContainerType3& weights,
            value_type eps,
            value_type nrmb_correction = 1.)
    {
        // Should this be another class?
        // if A does not change Lanczos iterations could be reused from apply function!?
#ifdef MPI_VERSION
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif //MPI
        dg::Timer t;
        t.tic();
        auto func = make_FuncEigen_Te1( [&](value_type x) {return op( x, 1.);});
        auto T = m_lanczos.tridiag( func, std::forward<MatrixType>(A),
                b, weights, eps, nrmb_correction,
                "universal", 1.0, 2);
        compute_vlcl_adjoint( op, std::forward<MatrixType>(A), diag, T, x, b,
                weights, m_lanczos.get_bnorm());

        t.toc();
        if( m_benchmark)
            DG_RANK0 std::cout << "# `"<<m_message<<"` solve with {"<<T.num_rows<<"} iterations took "<<t.diff()<<"s\n";
        return T.num_rows;
    }

    template< class BinaryOp, class ContainerType0, class MatrixType,
        class ContainerType1, class ContainerType2>
    void compute_vlcl( BinaryOp op, const ContainerType0& diag,
            MatrixType&& A,
            const TriDiagonal<thrust::host_vector<value_type>>& T,
            ContainerType1& x,
            const ContainerType2& b,
            value_type bnorm)
    {
        dg::blas1::copy(0., x);
        if( 0 == bnorm )
        {
            return;
        }
        unsigned iter = T.O.size();
        thrust::host_vector<value_type> evals = T.O, plus  = T.P;
        thrust::host_vector<value_type> work (2*iter-2);
        dg::SquareMatrix<value_type> EHt(iter);
        //Compute Eigendecomposition
        lapack::stev(LAPACK_COL_MAJOR, 'V', evals, plus, EHt.data(), work);
        dg::blas1::axpby(1./bnorm, b, 0.0, m_v); //m_v[1] = b/||b||
        dg::blas1::copy(0., m_vm);
        // compute c_1 v_1
        for ( unsigned k=0; k<iter; k++)
        {
            dg::blas1::evaluate( m_f, dg::equals(), op, diag, evals[k]);
            dg::blas1::pointwiseDot( bnorm*EHt(k, 0)*EHt(k,0), m_f, m_v, 1.,
                    x);
        }
        for ( unsigned i=0; i<iter-1; i++)
        {
            dg::blas2::symv( std::forward<MatrixType>(A), m_v, m_vp);
            dg::blas1::axpbypgz(
                    -T.M[i]/T.P[i], m_vm,
                    -T.O[i]/T.P[i], m_v,
                        1.0/T.P[i], m_vp);
            m_vm.swap( m_v);
            m_v.swap( m_vp);
            // compute c_l v_l
            for ( unsigned k=0; k<iter; k++)
            {
                dg::blas1::evaluate( m_f, dg::equals(), op, diag, evals[k]);
                dg::blas1::pointwiseDot( bnorm*EHt(k,0)*EHt(k,i+1), m_f, m_v,
                        1., x);
            }
        }
    }
    template< class BinaryOp, class MatrixType, class ContainerType0,
        class ContainerType1, class ContainerType2, class ContainerType3>
    void compute_vlcl_adjoint( BinaryOp op,
            MatrixType&& A,
            const ContainerType0& diag,
            const TriDiagonal<thrust::host_vector<value_type>>& T,
            ContainerType1& x,
            const ContainerType2& b,
            const ContainerType3& weights,
            value_type bnorm)
    {
        dg::blas1::copy(0., x);
        if( 0 == bnorm )
        {
            return;
        }
        unsigned iter = T.O.size();
        thrust::host_vector<value_type> evals = T.O, plus  = T.P;
        thrust::host_vector<value_type> work (2*iter-2);
        dg::SquareMatrix<value_type> EHt(iter);
        //Compute Eigendecomposition
        lapack::stev(LAPACK_COL_MAJOR, 'V', evals, plus, EHt.data(), work);
        dg::blas1::axpby(1./bnorm, b, 0.0, m_v); //m_v[1] = b/||b||
        dg::blas1::copy(0., m_vm);
        // compute alpha_i1
        dg::SquareMatrix<value_type> alpha(iter);
        for ( unsigned k=0; k<iter; k++)
        {
            dg::blas1::evaluate( m_f, dg::equals(), op, evals[k], diag);
            dg::blas1::pointwiseDot( m_f, m_v, m_f);
            alpha( k,0) = dg::blas2::dot( m_f, weights, b);
        }
        for ( unsigned i=0; i<iter-1; i++)
        {
            dg::blas2::symv( std::forward<MatrixType>(A), m_v, m_vp);
            dg::blas1::axpbypgz(
                    -T.M[i]/T.P[i], m_vm,
                    -T.O[i]/T.P[i], m_v,
                        1.0/T.P[i], m_vp);
            m_vm.swap( m_v);
            m_v.swap( m_vp);
            for ( unsigned k=0; k<iter; k++)
            {
                dg::blas1::evaluate( m_f, dg::equals(), op, evals[k], diag);
                dg::blas1::pointwiseDot( m_f, m_v, m_f);
                alpha( k,i+1) = dg::blas2::dot( m_f, weights, b);
            }
        }
        // Observation: With an exponential function the lines of alpha get extremely small (because exp(lambda) gets very small ... so maybe one can save a few scalar products
        // compute E_li E_ki alpha_ik v_l
        std::vector<double> cl( iter, 0.0);
        for( unsigned l=0; l<iter; l++)
            for( unsigned i=0; i<iter; i++)
                for( unsigned k=0; k<iter; k++)
                    cl[l] += EHt(k,i)*alpha(k,i)*EHt(k,l);
        // 3rd Lanczos iteration
        dg::blas1::axpby(1./bnorm, b, 0.0, m_v); //m_v[1] = b/||b||
        dg::blas1::copy(0., m_vm);
        dg::blas1::axpby( cl[0], m_v, 1., x);
        for ( unsigned i=0; i<iter-1; i++)
        {
            dg::blas2::symv( std::forward<MatrixType>(A), m_v, m_vp);
            dg::blas1::axpbypgz(
                    -T.M[i]/T.P[i], m_vm,
                    -T.O[i]/T.P[i], m_v,
                        1.0/T.P[i], m_vp);
            m_vm.swap( m_v);
            m_v.swap( m_vp);
            dg::blas1::axpby( cl[i+1], m_v, 1., x);
        }
    }
    UniversalLanczos<ContainerType>& lanczos() { return m_lanczos;}
    private:

    UniversalLanczos<ContainerType> m_lanczos;
    bool m_benchmark = true;
    std::string m_message = "ProductFunction";
    ContainerType  m_v, m_vp, m_vm, m_f;
};

}//namespace mat
}//namespace dg
