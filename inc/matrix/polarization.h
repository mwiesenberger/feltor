#pragma once
#include <memory>
#include "dg/algorithm.h"

#include "lanczos.h"
#include "matrixsqrt.h"
#include "matrixfunction.h"
#include "tensorelliptic.h"

namespace dg {
namespace mat {

/**
 * @brief Various arbitary wavelength polarization charge operators of delta-f
 *  (df) and full-f (ff)
 *
 * "df"
 * \f[
 * x = -\Delta \left(1+\alpha\Delta\right)^{-1} \phi  \f]  \f[
 * \f]
 * "ff"
 * \f[
 * x = \sqrt{1+\alpha\Delta}^{-1}\left(-\nabla \cdot \chi
 *  \nabla\right)\sqrt{1+\alpha\Delta}^{-1}\phi \f]
 * "ffO4"
 * \f[
 *   x = (1+\alpha\Delta)^{-1} \left(-\nabla \cdot \chi \nabla - \Delta \iota \Delta
 *   +  \nabla \cdot\nabla \cdot 2\iota \nabla \nabla \right)(1+\alpha\Delta)^{-1}
 *   \phi
 * \f]
 *
 * @ingroup matrixmatrixoperators
 *
 */
template <class Geometry, class Matrix, class Container>
class PolCharge
{
    public:
    using value_type = get_value_type<Container>;
    ///@brief empty object ( no memory allocation)
    PolCharge(){}
    /**
     * @brief Construct from Grid
     *
     * @param alpha alpha of the Helmholtz operator
     * @param eps_gamma epsilon (-vector) for the Helmholtz operator inversion
     *  or the sqrt Helmholtz operator inversion ("df" and "ffO4" need one
     *  entry per multigrid stage, i.e. 3, "ff" uses only \c eps_gamma[0])
     * @param g The Grid, boundary conditions are taken from here
     * @param dir Direction of the right first derivative in x and y
     *  (i.e. \c dg::forward, \c dg::backward or \c dg::centered),
     * @param jfactor (\f$ = \alpha \f$ ) scale jump terms (1 is a good value
     *  but in some cases 0.1 or 0.01 might be better)
     * @param mode arbitrary wavelength polarization charge mode ("df" / "ff" /
     *  "ffO4" are implemented, other values throw a \c dg::Error)
     * @param commute false if Helmholtz operators (or their square root) are
     *  outside the elliptic or tensorelliptic operator and true otherwise
     *  (true is only implemented for "df", \c symv throws a \c dg::Error otherwise)
     * @note For \c alpha=0 the Helmholtz operators are the identity and only
     *  the elliptic or tensor elliptic operator is applied
     */
    PolCharge(value_type alpha, std::vector<value_type> eps_gamma,
            const Geometry& g, direction dir = forward, value_type jfactor=1.,
            std::string mode = "df", bool commute = false):
        PolCharge( alpha, eps_gamma, g, g.bcx(), g.bcy(), dir, jfactor, mode,
                commute)
    { }
    /**
     * @brief Construct from boundary conditions
     *
     * @param alpha alpha of the Helmholtz operator
     * @param eps_gamma epsilon (-vector) for the Helmholtz operator inversion
     *  or the sqrt Helmholtz operator inversion ("df" and "ffO4" need one
     *  entry per multigrid stage, i.e. 3, "ff" uses only \c eps_gamma[0])
     * @param g The Grid, boundary conditions are taken from here
     * @param bcx boundary condition in x
     * @param bcy boundary contition in y
     * @param dir Direction of the right first derivative in x and y
     *  (i.e. \c dg::forward, \c dg::backward or \c dg::centered),
     * @param jfactor (\f$ = \alpha \f$ ) scale jump terms (1 is a good value
     *  but in some cases 0.1 or 0.01 might be better)
     * @param mode arbitrary wavelength polarization charge mode ("df" / "ff" /
     *  "ffO4" (O as in Order!) are implemented, other values throw a \c dg::Error)
     * @param commute false if Helmholtz operators (or their square root) are
     *  outside the elliptic or tensorelliptic operator and true otherwise
     *  (true is only implemented for "df", \c symv throws a \c dg::Error otherwise)
     * @note For \c alpha=0 the Helmholtz operators are the identity and only
     *  the elliptic or tensor elliptic operator is applied
    */
    PolCharge( value_type alpha, std::vector<value_type> eps_gamma,
            const Geometry& g, bc bcx, bc bcy, direction dir = forward,
            value_type jfactor=1., std::string mode = "df", bool commute = false)
    {
        m_alpha = alpha;
        m_eps_gamma = eps_gamma;
        m_mode = mode;
        m_commute = commute;
        if( m_mode != "df" && m_mode != "ff" && m_mode != "ffO4")
            throw dg::Error( dg::Message(_ping_)<<"PolCharge mode "<<m_mode<<" not implemented!");
        m_temp2 = dg::evaluate(dg::zero, g);
        m_temp =  m_temp2;
        m_temp2_ex.set_max(1, m_temp2);
        m_temp_ex.set_max(1, m_temp);
        if (m_mode == "df")
        {
            m_ell.construct(g, bcx, bcy, dir, jfactor );
            m_multi_g.construct(g, 3);
            for( unsigned u=0; u<3; u++)
            {
                m_multi_gamma.push_back( {m_alpha, {m_multi_g.grid(u), bcx, bcy,
                        dir, jfactor}});
            }
        }
        if (m_mode == "ff")
        {
            m_ell.construct(g, bcx, bcy, dir, jfactor );
            // m_inv_sqrt stores the Helmholtz operator by reference, so keep
            // it on the heap where copies of *this can share it
            // (for alpha = 0 the operator is the identity and not needed)
            if( m_alpha != 0)
            {
                m_ff_gamma = std::make_shared<dg::Helmholtz<Geometry, Matrix,
                    Container>>( m_alpha, dg::Elliptic<Geometry, Matrix,
                    Container>{g, bcx, bcy, dir, jfactor});
                m_inv_sqrt.construct( *m_ff_gamma, -1, g, m_eps_gamma[0]);
            }
        }
        if (m_mode == "ffO4")
        {
            m_tensorell.construct(g, bcx, bcy, dir, jfactor);
            m_multi_g.construct(g, 3);
            for( unsigned u=0; u<3; u++)
            {
                m_multi_gamma.push_back({ m_alpha, {m_multi_g.grid(u), bcx, bcy,
                        dir, jfactor}});
            }
        }
    }

    ///@copydoc hide_construct
    template<class ...Params>
    void construct( Params&& ...ps)
    {
        //construct and swap
        *this = PolCharge( std::forward<Params>( ps)...);
    }
    /**
     * @brief Change \f$\chi\f$ in the elliptic or tensor elliptic operator
     *
     * @param sigma The new scalar part in \f$\chi\f$ (all elements must be >0)
     * @tparam ContainerType0 must be usable with \c Container in \ref dispatch
     */
    template<class ContainerType0>
    void set_chi( const ContainerType0& sigma)
    {
        if (m_mode == "ff")
            m_ell.set_chi(sigma);
        if (m_mode == "ffO4")
        {
            m_tensorell.set_chi(sigma);
        }
    }
    /**
     * @brief Change \f$\chi\f$ in the elliptic or tensor elliptic operator
     *
     * @param tau The new tensor part in \f$\chi\f$ (all elements must be >0)
     * @tparam ContainerType0 must be usable in \c dg::assign to \c Container
     */
    template<class ContainerType0>
    void set_chi( const SparseTensor<ContainerType0>& tau)
    {
        if (m_mode == "ff")
            m_ell.set_chi(tau);
        if (m_mode == "ffO4")
        {
            m_tensorell.set_chi(tau);
        }
    }
    /**
     * @brief Change \f$\iota\f$ in the tensor elliptic operator
     *
     * @param sigma The new \f$\iota\f$ (all elements must be >0; only used in "ffO4")
     * @tparam ContainerType0 must be usable with \c Container in \ref dispatch
     */
    template<class ContainerType0>
    void set_iota( const ContainerType0& sigma)
    {
        if (m_mode == "ffO4")
        {
            m_tensorell.set_iota(sigma);
        }
    }
    /**
     * @brief Set the commute
     * @param commute Either true or false.
     */
    void set_commute( bool commute) {m_commute = commute;}
    /**
     * @brief Get the current state of commute
     * @return Either true or false.
     */
    bool get_commute() const {return m_commute;}
    /**
     * @brief Return the vector making the matrix symmetric
     *
     * i.e. the volume form
     * @return volume form including weights
     */
    const Container& weights()const {
        if (m_mode == "ffO4")
             return  m_tensorell.weights();
        else return  m_ell.weights();
    }
    /**
     * @brief Return the default preconditioner to use in conjugate gradient
     *
     * Returns the preconditioner of the elliptic (or, in "ffO4", tensor
     * elliptic) operator
     * @return preconditioner
     */
    const Container& precond()const {
        if (m_mode == "ffO4")
            return m_tensorell.precond();
        else return m_ell.precond();
    }

    /**
     * @brief Compute elliptic term and store in output
     *
     * i.e. \c y=M*x
     * @param x left-hand-side
     * @param y result
     * @tparam ContainerTypes must be usable with \c Container in \ref dispatch
     */
    template<class ContainerType0, class ContainerType1>
    void operator()( const ContainerType0& x, ContainerType1& y){
        symv( 1, x, 0, y);
    }

    /**
     * @brief Compute elliptic term and store in output
     *
     * i.e. \c y=M*x
     * @param x left-hand-side
     * @param y result
     * @tparam ContainerTypes must be usable with \c Container in \ref dispatch
     */
    template<class ContainerType0, class ContainerType1>
    void symv( const ContainerType0& x, ContainerType1& y){
        symv( 1, x, 0, y);
    }
    /**
     * @brief Compute elliptic term and add to output
     *
     * i.e. \c y=alpha*M*x+beta*y
     * @param alpha a scalar
     * @param x left-hand-side
     * @param beta a scalar
     * @param y result
     * @tparam ContainerTypes must be usable with \c Container in \ref dispatch
     */
    template<class ContainerType0, class ContainerType1>
    void symv( value_type alpha, const ContainerType0& x, value_type beta, ContainerType1& y)
    {

        if ( m_commute && m_mode != "df")
            throw dg::Error( dg::Message(_ping_)<<"PolCharge mode "<<m_mode
                    <<" with commute = true is not implemented!");
        if (m_alpha == 0) // the Helmholtz operators are the identity
        {
            if (m_mode == "ffO4")
                m_tensorell.symv(alpha, x, beta, y);
            else
                m_ell.symv(alpha, x, beta, y);
            return;
        }
        if (m_mode == "df")
        {
            if (m_commute == false)
            {
                m_temp2_ex.extrapolate(m_temp2);
                std::vector<unsigned> number = m_multi_g.solve(
                        m_multi_gamma, m_temp2, x, m_eps_gamma);
                m_temp2_ex.update(m_temp2);

                m_ell.symv(alpha, m_temp2, beta, y);
            }
            else
            {
                m_ell.symv(1.0, x, 0.0, m_temp);

                m_temp2_ex.extrapolate(m_temp2);
                std::vector<unsigned> number = m_multi_g.solve(
                        m_multi_gamma, m_temp2, m_temp, m_eps_gamma);
                m_temp2_ex.update(m_temp2);

                dg::blas1::axpby(alpha, m_temp2, beta, y);

            }

        }
        if (m_mode == "ff" ) //assuming constant FLR effects
        {
            if (m_commute == false)
            {
                //unsigned number = 0 ;
                dg::apply( m_inv_sqrt, x, m_temp2);
                //std::cout << "#number of sqrt iterations: " << number << " "<<m_eps_gamma[0]<< std::endl;

                m_ell.symv(1.0, m_temp2, 0.0, m_temp);

                dg::apply( m_inv_sqrt, m_temp, m_temp2);
                //std::cout << "#number of sqrt iterations: " << number << std::endl;
                //number++;//avoid compiler warning

                dg::blas1::axpby(alpha, m_temp2, beta, y);
            }
        }
        if (m_mode == "ffO4")
        {
            if (m_commute == false)
            {
                m_temp2_ex.extrapolate(m_temp2);
                std::vector<unsigned> number = m_multi_g.solve(
                        m_multi_gamma, m_temp2, x, m_eps_gamma);
                m_temp2_ex.update(m_temp2);

                m_tensorell.symv(1.0, m_temp2, 0.0, m_temp);

                m_temp_ex.extrapolate(m_temp2);
                number = m_multi_g.solve( m_multi_gamma, m_temp2,
                        m_temp, m_eps_gamma);
                m_temp_ex.update(m_temp2);

                dg::blas1::axpby(alpha, m_temp2, beta, y);
            }
        }
    }

    private:
    dg::Elliptic<Geometry,  Matrix, Container> m_ell;
    dg::mat::TensorElliptic<Geometry,  Matrix, Container> m_tensorell;

    std::vector< dg::Helmholtz<Geometry,  Matrix, Container> > m_multi_gamma;
    std::shared_ptr< dg::Helmholtz<Geometry,  Matrix, Container> > m_ff_gamma;
    dg::MultigridCG2d<Geometry, Matrix, Container> m_multi_g;
    dg::mat::MatrixSqrt<Container> m_inv_sqrt;
    Container m_temp, m_temp2;
    value_type  m_alpha = 0;
    std::vector<value_type> m_eps_gamma;
    std::string m_mode;
    dg::Extrapolation<Container> m_temp2_ex, m_temp_ex;
    bool m_commute = false;


};

}  //namespace mat

///@cond
template< class G, class M, class V>
struct TensorTraits< mat::PolCharge<G, M, V> >
{
    using value_type      = get_value_type<V>;
    using tensor_category = SelfMadeMatrixTag;
};
///@endcond
}  //namespace dg
