#pragma once

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "parameters.h"

namespace thermal
{
template< class Geometry, class IMatrix, class Matrix, class Container >
struct Collisions
{
    Collisions( const Geometry& g, thermal::Parameters p,
        dg::geo::TokamakMagneticField mag, dg::file::WrappedJsonValue js):m_p(p), m_js(js)
    {
        dg::assign( dg::evaluate( dg::zero, g), m_temp0 );
        m_temp1 = m_temp0;
        dg::assign(  dg::pullback(dg::geo::Bmodule(mag), g), m_B2);
        dg::blas1::pointwiseDot( m_B2, m_B2, m_B2);
        m_tperpST.resize( p.num_species);
        m_tparaST.resize( p.num_species);
        m_u.resize( p.num_species);
        m_lapperpN.construct ( g, p.bcxN, p.bcyN,  p.diff_dir);
        m_R0 = mag.R0();
    }

    void gather_quantities( unsigned s, const std::array<Container, 6>& q, const std::array<Container,6>& qST)
    {
        dg::blas1::copy( qST[1], m_tperpST[s]);
        dg::blas1::copy( qST[2], m_tparaST[s]);
        dg::blas1::copy( q[3], m_u[s]);
    }

    // After calling gather_quantities for all s:
    const std::vector<Container>& get_velocity() const{ return m_u;}
    const std::vector<Container>& get_tperpST() const{ return m_tperpST;}
    const std::vector<Container>& get_tparaST() const{ return m_tparaST;}

    void add_collisions(
        const std::vector<Container>& densityST,
        const Container& aparST,
        const std::array<std::vector<Container>,6>& y,
        std::array<std::vector<Container>,6>& yp);
    private:
    const thermal::Parameters m_p;
    const dg::file::WrappedJsonValue m_js;
    std::vector<Container> m_tperpST, m_tparaST, m_u; // gather through species loop
    dg::Elliptic2d< Geometry, Matrix, Container> m_lapperpN;
    Container m_temp0, m_temp1, m_B2;
    double m_R0;
};
template<class Grid, class IMatrix, class Matrix, class Container>
void Collisions<Grid, IMatrix, Matrix, Container>::add_collisions(
        const std::vector<Container>& densityST,
        const Container& aparST,
        const std::array<std::vector<Container>,6>& y,
        std::array<std::vector<Container>,6>& yp)
{
    double nu_ref = m_p.nu_ref;
    for( unsigned s=0; s<m_p.num_species; s++)
    {
        double zs = m_p.z[s], mus = m_p.mu[s];
        // spperp Coulomb
        dg::blas1::copy( 0., m_temp0);
        for( unsigned k=0; k<m_p.num_species; k++)
            if( s != k)
            {
                double zk = m_p.z[k];
                double muk = m_p.mu[k];
                dg::blas1::subroutine( [nu_ref, zs, zk, mus, muk] DG_DEVICE(
                    double& spperp, double ns, double nk, double ps, double pk)
                {
                    double ts = ps/ns, tk = pk/nk;
                    spperp += -2*nu_ref*ns*nk*zs*zs*zk*zk/mus/muk/(
                        ts/mus + tk/muk)/sqrt( ts/mus + tk/muk)*(ts-tk);
                }, m_temp0, y[0][s], y[0][k], y[1][s], y[1][k]);
            }
        dg::blas1::axpby( 1., m_temp0, 1., yp[1][s]);
        // spperp Lorentz
        double pis = m_p.pi[s];
        dg::blas1::subroutine( [nu_ref, zs, mus, pis] DG_DEVICE(
            double & spperp, double ns, double pperps, double pparas)
        {
            double ts = pperps/ns;
            double nu_ss = nu_ref*ns*zs*zs*zs*zs/sqrt(2*mus*ts)/ts;
            spperp += nu_ss/pis/3*(pparas-pperps);
        }, yp[1][s], y[0][s], y[1][s], y[2][s]);


        // sn
        dg::blas1::pointwiseDivide( m_p.mu[s]/2./m_p.z[s]/m_p.z[s], m_temp0, m_B2, 0., m_temp0);
        dg::blas2::symv( m_lapperpN, m_temp0, m_temp1); // should have 0 BC due to (tperp - tpara)
        dg::blas1::axpby( 1., m_temp1, 1., yp[0][s]);

        // sppara Coulomb
        dg::blas1::pointwiseDot( m_p.mu[s], m_u[s], m_u[s], m_temp1, 0., m_temp0); // mu U^2 S_N
        for( unsigned k=0; k<m_p.num_species; k++)
            if( s != k)
            {
                double zs = m_p.z[s], zk = m_p.z[k];
                double mus = m_p.mu[s], muk = m_p.mu[k];
                dg::blas1::subroutine( [nu_ref, zs, zk, mus, muk] DG_DEVICE(
                    double& sppara, double ns, double nk, double ps, double pk, double us, double uk)
                {
                    double ts = ps/ns, tk = pk/nk;
                    sppara += -2*nu_ref*ns*nk*zs*zs*zk*zk/mus/muk/(
                        ts/mus + tk/muk)/sqrt( ts/mus + tk/muk)*((ts-tk) - 0.51*muk*(us-uk)*(us-uk));
                }, m_temp0, y[0][s], y[0][k], y[2][s], y[2][k], m_u[s], m_u[k]);
            }
        dg::blas1::axpby( 1., m_temp0, 1., yp[2][s]);
        // sppara Lorentz
        dg::blas1::subroutine( [nu_ref, zs, mus,pis] DG_DEVICE(
            double & sppara, double ns, double pperps, double pparas)
        {
            double ts = pparas/ns;
            double nu_ss = nu_ref*ns*zs*zs*zs*zs/sqrt(2*mus*ts)/ts;
            sppara += -2.*nu_ss/pis/3*(pparas-pperps);
        }, yp[2][s], y[0][s], y[1][s], y[2][s]);


        // SU Coulomb
        // spperpST Coulomb
        dg::blas1::copy( 0., m_temp0);
        for( unsigned k=0; k<m_p.num_species; k++)
        if( s != k)
        {
            double zs = m_p.z[s], zk = m_p.z[k];
            double mus = m_p.mu[s], muk = m_p.mu[k];
            dg::blas1::subroutine( [nu_ref, zs, zk, mus, muk] DG_DEVICE(
                double& spperp, double ns, double nk, double ts, double tk)
            {
                spperp += -2*nu_ref*ns*nk*zs*zs*zk*zk/mus/muk/(
                    ts/mus + tk/muk)/sqrt( ts/mus + tk/muk)*(ts-tk);
            }, m_temp0, densityST[s], densityST[k], m_tperpST[s], m_tperpST[k]);
        }
        dg::blas1::pointwiseDivide( m_p.mu[s]/2./m_p.z[s]/m_p.z[s], m_temp0, m_B2, 0., m_temp0);
        dg::blas2::symv( m_lapperpN, m_temp0, m_temp1); // S_NST
        // SU ST
        dg::blas1::axpby( 1., y[3][s], -m_p.z[s]/m_p.mu[s], aparST, m_temp0);
        dg::blas1::pointwiseDivide( m_temp0, y[0][s], m_temp0); // U/N
        dg::blas1::pointwiseDot( -1., m_temp0, m_temp1, 0., m_temp0); // -U/N S_NST
        for( unsigned k=0; k<m_p.num_species; k++)
        if( s != k)
        {
            double zs = m_p.z[s], zk = m_p.z[k];
            double mus = m_p.mu[s], muk = m_p.mu[k];
            dg::blas1::subroutine( [nu_ref, zs, zk, mus, muk] DG_DEVICE(
                double& su, double nk, double ts, double tk, double ws, double wk, double apar)
            {
                double us = ws - zs/mus*apar, uk = wk - zk/muk*apar;
                su += -0.51*nu_ref*nk*zs*zs*zk*zk*(mus+muk)/mus/mus/muk/(
                    ts/mus + tk/muk)/sqrt( ts/mus + tk/muk)*(us-uk);
            }, m_temp0, densityST[k], m_tparaST[s], m_tparaST[k], y[3][s], y[3][k], aparST);
        }

        // SQperp Lorentz
        double kappas = m_p.kappa[s];
        double qlandau = m_p.qlandau;
        double R0 = m_R0;
        dg::blas1::subroutine( [nu_ref, zs, mus,kappas, qlandau, R0] DG_DEVICE(
            double & sqperp, double ns, double ts, double qperp, double qpara)
        {
            double nu_ss = nu_ref*ns*zs*zs*zs*zs/sqrt(2*mus*ts)/ts;
            sqperp += - 5/2./kappas*nu_ss*( qperp - 1.28*(0.5*qpara - 1.5*qperp))
                      - 1./qlandau/R0*sqrt( ts/mus)*qperp;
        }, yp[4][s], densityST[s], m_tperpST[s], y[4][s], y[5][s]);

        // SQpara Lorentz
        dg::blas1::subroutine( [nu_ref, zs, mus,kappas, qlandau, R0] DG_DEVICE(
            double & sqpara, double ns, double ts, double qperp, double qpara)
        {
            double nu_ss = nu_ref*ns*zs*zs*zs*zs/sqrt(2*mus*ts)/ts;
            sqpara += - 5/2./kappas*nu_ss*( qperp - 1.28*(0.5*qpara - 1.5*qperp))
                      - 1./qlandau/R0*sqrt( ts/mus)*qperp;
        }, yp[5][s], densityST[s], m_tparaST[s], y[4][s], y[5][s]);

    }

}

}
