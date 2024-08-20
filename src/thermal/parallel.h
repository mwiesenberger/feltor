#pragma once

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "parameters.h"

#define FELTORPARALLEL 1

namespace thermal
{
template< class Geometry, class IMatrix, class Matrix, class Container >
class ParallelDynamics
{
    ParallelDynamics( const Geometry&, thermal::Parameters,
        dg::geo::TokamakMagneticField, dg::file::WrappedJsonValue);
    const dg::geo::Fieldaligned<Geometry, IMatrix, Container>& fieldaligned() const
    {
        return m_fa;
    }
    // call before computing AparST
    void udpate_staggered_density( const std::vector<Container>& density)
    {
        for( unsigned s=0; s<m_p.num_species; s++)
        {
            //density
            m_faST( dg::geo::zeroMinus, density[s], m_minusSTN[s]);
            m_faST( dg::geo::einsPlus,  density[s], m_STN[s]);
            update_parallel_bc_1st( m_minusSTN[s], m_STN[s],
                    m_p.bcxN, m_p.bcxN == dg::DIR ? m_p.nbc[s] : 0.);
            dg::blas1::axpby( 0.5, m_minusSTN[s], 0.5, m_STN[s], m_STN[s]);
        }
    }
    const std::vector<Container>& get_staggered_density() const{
        return m_STN;
    }
    // only call once
    void update_apar( const Container& aparST)
    {
        m_faST( dg::geo::einsMinus, aparST, m_minusST[0]);
        m_faST( dg::geo::zeroPlus,  aparST, m_plusST[0]);
        update_parallel_bc_1st( m_minusST[0], m_plusST[0], m_p.bcxA, 0.);
        dg::blas1::axpby( 0.5, m_minusST[0], 0.5, m_plusST[0], m_apar);
    }
    const Container& get_apar() const{ return m_apar;}

    // call once per species
    void update_quantities( unsigned s, const Container& aparST, const std::array<std::vector<Container>,6>& y);
    // N, Tperp, Tpara, U, Uperp, Upara
    const std::array<Container,6>& get_q( ) const{return m_q;}
    const std::array<Container,6>& get_qST( ) const{return m_qST;}

    private:
    dg::geo::Fieldaligned<Geometry, IMatrix, Container> m_fa, m_faST;

    std::vector<Container> m_minusSTN, m_STN;
    // m_plusSTN = 2*m_STN - m_minusSTN

    // Everything is:
    // N, Tperp, Tpara, U, Uperp, Upara
    std::array<Container,6> m_minus, m_zero, m_plus;
    std::array<Container,6> m_minusST, m_plusST, m_q, m_qST;

    std::array<Container,4> m_minusSTGamma, m_plusSTGamma;
    Container m_temp, m_apar, m_gradParGamma2;

    const thermal::Parameters m_p;
    const dg::file::WrappedJsonValue m_js;
    bool m_reversed_field = false;
};

template<class Grid, class IMatrix, class Matrix, class Container>
ParallelDynamics<Grid, IMatrix, Matrix, Container>::ParallelDynamics( const Grid& g,
    thermal::Parameters p, dg::geo::TokamakMagneticField mag,
    dg::file::WrappedJsonValue js
    ): m_p(p), m_js(js)
{
    dg::assign( dg::evaluate( dg::zero, g), m_temp );
    m_minusSTN.resize( m_p.num_species);
    std::fill( m_minusSTN.begin(), m_minusSTN.end(), m_temp);
    m_STN = m_minusSTN;
    std::fill( m_minus.begin(), m_minus.end(), m_temp);
    m_zero = m_plus = m_minus;
    std::fill( m_minusST.begin(), m_minusST.end(), m_temp);
    std::fill( m_plusST.begin(), m_plusST.end(), m_temp);
    m_q = m_qST = m_plusST;

    std::fill( m_plusSTGamma.begin(), m_plusSTGamma.end(), m_temp);
    m_apar = m_gradParGamma2 = m_temp;

    m_reversed_field = false;
    if( mag.ipol()( g.x0(), g.y0()) < 0)
        m_reversed_field = true;
    //in DS we take the true bhat
    auto bhat = dg::geo::createBHat( mag);
    // do not construct FCI if we just want to calibrate
    if( !p.calibrate )
    {
        m_fa.construct( bhat, g, dg::NEU, dg::NEU, dg::geo::NoLimiter(),
            p.rk4eps, p.mx, p.my, 2.*M_PI/(double)p.Nz, p.interpolation_method);
        m_faST.construct( bhat, g, dg::NEU, dg::NEU, dg::geo::NoLimiter(),
            p.rk4eps, p.mx, p.my, 2.*M_PI/(double)p.Nz/2., p.interpolation_method );
    }
}
template<class Grid, class IMatrix, class Matrix, class Container>
void ParallelDynamics<Grid, IMatrix, Matrix, Container>::update_quantities(
        unsigned s,
        const Container& aparST,
        const std::array<std::vector<Container>,6>& y)
{
#ifdef MPI_VERSION
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
    // 1st transform all densities
    std::vector<Container>&  density    = y[0];
    std::vector<Container>&  pperp      = y[1];
    std::vector<Container>&  ppara      = y[2];
    // N, Tperp, Tpara
    dg::blas1::copy( density[s], m_q[0]);
    dg::blas1::pointwiseDivide( pperp[s], density[s], m_q[1]);
    dg::blas1::pointwiseDivide( ppara[s], density[s], m_q[2]);
    std::array<double, 3> nbc = {m_p.nbc[s], m_p.tbc, m_p.tbc};
    for( unsigned u=0; u<3; u++)
    {
        m_fa( dg::geo::einsMinus, m_q[u], m_minus[u]);
        m_fa( dg::geo::zeroForw,  m_q[u], m_zero[u]);
        m_fa( dg::geo::einsPlus,  m_q[u], m_plus[u]);
        update_parallel_bc_2nd( m_fa, m_minus[u], m_zero[u], m_plus[u],
                m_p.bcxN, m_p.bcxN == dg::DIR ? nbc[u] : 0.);
        if( u!= 0)
        {
            m_faST( dg::geo::zeroMinus, m_q[s], m_minusST[u]);
            m_faST( dg::geo::einsPlus,  m_q[s], m_plusST[u]);
            update_parallel_bc_1st( m_minusST[u], m_plusST[u],
                m_p.bcxN, m_p.bcxN == dg::DIR ? nbc[u] : 0.);
            dg::blas1::axpby( 0.5, m_minusST[u], 0.5, m_plusST[u], m_qST[u]);
        }
        else
        {
            dg::blas1::copy( m_minusSTN[s], m_minusST[0]);
            dg::blas1::axpby( 2., m_STN[s], -1., m_minusSTN[s], m_plusST[0]);
            dg::blas1::copy( m_STN[s], m_qST[0]);
        }
    }
    // 2nd transform all velocities
    std::vector<Container>&  wST        = y[3];
    std::vector<Container>&  qperpST    = y[4];
    std::vector<Container>&  qparaST    = y[5];
    // UST, UperpST, UparaST
    dg::blas1::axpby( 1., wST[s], -m_p.z[s]/m_p.mu[s], aparST, m_qST[3]);
    dg::blas1::pointwiseDot( m_qST[0], m_qST[1], m_temp); // pperpST
    dg::blas1::pointwiseDivide( qperpST[s], m_temp, m_qST[4]);
    dg::blas1::pointwiseDot( m_qST[0], m_qST[2], m_temp); // pparaST
    dg::blas1::pointwiseDivide( qparaST[s], m_temp, m_qST[5]);
    for( unsigned u=3; u<6; u++)
    {
        m_fa( dg::geo::einsMinus, m_qST[u], m_minus[u]);
        m_fa( dg::geo::zeroForw,  m_qST[u], m_zero[u]);
        m_fa( dg::geo::einsPlus,  m_qST[u], m_plus[u]);
        update_parallel_bc_2nd( m_fa, m_minus[u], m_zero[u],
            m_plus[u], m_p.bcxU, 0.);
        m_faST( dg::geo::einsMinus, m_qST[u], m_minusST[u]);
        m_faST( dg::geo::zeroPlus,  m_qST[u],  m_plusST[u]);
        update_parallel_bc_1st( m_minusST[u], m_plusST[u], m_p.bcxU, 0.);
        dg::blas1::axpby( 0.5, m_minusST[u], 0.5, m_plusST[u], m_q[u]);
    }
    // Now we have q, qST, qminus, qzero, qplus, qminusST, qplusST
}

}//namespace thermal
