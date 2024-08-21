#pragma once

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "parameters.h"

#define FELTORPERP 1

namespace thermal
{


template< class Geometry, class IMatrix, class Matrix, class Container >
class PerpDynamics
{
    PerpDynamics( const Geometry&, thermal::Parameters,
        dg::geo::TokamakMagneticField, dg::file::WrappedJsonValue);

    void update_density_derivatives(
        unsigned s,
        const Container& apar,
        const std::array<std::vector<Container>,6> >& y,
        const std::array<Container,6>& q
    );
    void add_perp_density_dynamics(
        unsigned s,
        const Container& apar,
        const std::array<std::vector<Container>,6> >& y,
        const std::array<Container,6>& q,
        std::array<std::vector<Container>,6>& yp
    );
    void add_perp_velocity_dynamics(
        unsigned s,
        const Container& apar,
        const std::array<std::vector<Container>,6> >& y,
        const std::array<Container,6>& q,
        std::array<std::vector<Container>,6>& yp
    );
    private:
    //these should be considered const
    std::array<Container,2> m_curvNabla, m_curvKappa, m_gradLnB;
    Container m_divCurvKappa, m_b_2, m_divb; //m_bphiB is bhat/ sqrt(g) / B (covariant component)

    Matrix m_dxF_N, m_dxB_N, m_dxF_U, m_dxB_U, m_dx_N, m_dx_U, m_dx_P, m_dx_A;
    Matrix m_dyF_N, m_dyB_N, m_dyF_U, m_dyB_U, m_dy_N, m_dy_U, m_dy_P, m_dy_A;

    Container m_temp0, m_temp1;
    std::array<Container,3> m_dxF, m_dxB, m_dyF, m_dyB;
    std::array<Container,9> m_dx, m_dy;
    std::map<std::string, Container> m_v;

    const thermal::Parameters m_p;
    const dg::file::WrappedJsonValue m_js;
};

template<class Grid, class IMatrix, class Matrix, class Container>
PerpDynamics<Grid, IMatrix, Matrix, Container>::PerpDynamics( const Grid& g,
    thermal::Parameters p, dg::geo::TokamakMagneticField mag,
    dg::file::WrappedJsonValue js
    ):
    m_dxF_N( dg::create::dx( g, p.bcxN, dg::forward) ),
    m_dxB_N( dg::create::dx( g, p.bcxN, dg::backward) ),
    m_dxF_U( dg::create::dx( g, p.bcxU, dg::forward) ),
    m_dxB_U( dg::create::dx( g, p.bcxU, dg::backward) ),
    m_dx_N(  dg::create::dx( g, p.bcxN, dg::centered) ),
    m_dx_U(  dg::create::dx( g, p.bcxU, dg::centered) ),
    m_dx_P(  dg::create::dx( g, p.bcxP, p.pol_dir) ),
    m_dx_A(  dg::create::dx( g, p.bcxA, p.pol_dir) ),
    m_dyF_N( dg::create::dy( g, p.bcyN, dg::forward) ),
    m_dyB_N( dg::create::dy( g, p.bcyN, dg::backward) ),
    m_dyF_U( dg::create::dy( g, p.bcyU, dg::forward) ),
    m_dyB_U( dg::create::dy( g, p.bcyU, dg::backward) ),
    m_dy_N(  dg::create::dy( g, p.bcxN, dg::centered) ),
    m_dy_U(  dg::create::dy( g, p.bcxU, dg::centered) ),
    m_dy_P(  dg::create::dy( g, p.bcyP, p.pol_dir) ),
    m_dy_A(  dg::create::dy( g, p.bcyA, p.pol_dir) ),
    m_p(p), m_js(js)
{
    //--------------------------Construct-------------------------//
    //due to the various approximations bhat and mag not always correspond
    dg::geo::CylindricalVectorLvl0 curvNabla, curvKappa;
    bool reversed_field = false;
    if( mag.ipol()( g.x0(), g.y0()) < 0)
        reversed_field = true;
    if( p.curvmode == "true" )
        throw std::runtime_error( "curvmode : true is not possible in thermal code!");
    else if( p.curvmode == "low beta")
    {
        if( reversed_field)
            curvNabla = curvKappa = dg::geo::createCurvatureNablaB(mag, -1);
        else
            curvNabla = curvKappa = dg::geo::createCurvatureNablaB(mag, +1);
        dg::assign( dg::evaluate(dg::zero, g), m_divCurvKappa);
    }
    else if( p.curvmode == "toroidal")
    {
        if( reversed_field)
        {
            curvNabla = dg::geo::createCurvatureNablaB(mag, -1);
            curvKappa = dg::geo::createCurvatureKappa(mag, -1);
            dg::assign(  dg::pullback(dg::geo::DivCurvatureKappa(mag, -1), g),
                m_divCurvKappa);
        }
        else
        {
            curvNabla = dg::geo::createCurvatureNablaB(mag, +1);
            curvKappa = dg::geo::createCurvatureKappa(mag, +1);
            dg::assign(  dg::pullback(dg::geo::DivCurvatureKappa(mag, +1), g),
                m_divCurvKappa);
        }
    }
    else
        throw std::runtime_error( "Warning! curvmode value '"+p.curvmode+"' not recognized!! I don't know what to do! I exit!\n");
    dg::pushForward(curvNabla.x(), curvNabla.y(), curvNabla.z(),
        m_curvNabla[0], m_curvNabla[1], m_temp0, g);
    dg::pushForward(curvKappa.x(), curvKappa.y(), curvKappa.z(),
        m_curvKappa[0], m_curvKappa[1], m_temp0, g);
    dg::assign(  dg::pullback(dg::geo::Divb(mag), g), m_divb);
    // in PerpDynamics we take EPhi
    auto bhat = dg::geo::createEPhi(+1);
    if( p.curvmode == "true")
        throw std::runtime_error( "curvmode : true is not possible in thermal code!");
    else if( reversed_field)
        bhat = dg::geo::createEPhi(-1);
    dg::pushForward(bhat.x(), bhat.y(), bhat.z(), m_temp0, m_temp1, m_b_2, g);
    // make bhat covariant:
    dg::tensor::inv_multiply3d( g.metric(), m_temp0, m_temp1, m_b_2,
                                            m_temp0, m_temp1, m_b_2);
    // Grad Ln B covariant components
    m_gradLnB[0] = dg::pullback( dg::geo::BR(mag), g);
    m_gradLnB[1] = dg::pullback( dg::geo::BZ(mag), g);
    m_temp0 = dg::pullback(dg::geo::InvB(mag), g);
    dg::blas1::pointwiseDot( m_gradLnB[0], m_temp0, m_gradLnB[0]);
    dg::blas1::pointwiseDot( m_gradLnB[1], m_temp0, m_gradLnB[1]);
    dg::blas1::pointwiseDot( m_b_2, m_temp0, m_bphiB); // bphi / B
    m_temp0 = dg::tensor::volume( metric);
    dg::blas1::pointwiseDivide( m_b_2, m_temp0, m_bbphiB); //bphiB/detg/B

    // allocate resources for internal vectors
    std::fill( m_dxF.begin(), m_dxF.end(), m_temp0);
    m_dxB = m_dyF = m_dyB = m_dxF;
    std::fill( m_dx.begin(), m_dx.end(), m_temp0);
    std::fill( m_dy.begin(), m_dy.end(), m_temp0);
}

template<class Grid, class IMatrix, class Matrix, class Container>
void PerpDynamics<Grid, IMatrix, Matrix, Container>::update_density_derivatives(
    unsigned s,
    const Container& apar,
    const std::array<Container,4>& psi,
    const std::array<std::vector<Container>,6> >& y,
    const std::array<Container,6>& q
)
{
    const Container& to_advect[3] = {y[0][s], y[1][s], y[2][s]};
    double nbc[3] = {m_p.nbc[s], m_p.nbc[s]*m_p.tbc, m_p.nbc[s]*m_p.tbc};
    for( unsigned u=0; u<3; u++)
    {
        dg::blas1::transform( to_advect[u], m_temp0, dg::PLUS<double>(-nbc[u]));
        dg::blas2::symv( m_dxF_N, m_temp0, m_dxF[u]);
        dg::blas2::symv( m_dxB_N, m_temp0, m_dxB[u]);
        dg::blas2::symv( m_dyF_N, m_temp0, m_dyF[u]);
        dg::blas2::symv( m_dyB_N, m_temp0, m_dyB[u]);
    }
    // const Container& to_derive[9] = {apar, q[1], q[2], q[3], q[4], q[5], psi[0], psi[1], psi[2]};
    dg::blas2::symv( m_dx_A, apar, m_dx[0]);
    dg::blas2::symv( m_dy_A, apar, m_dy[0]);
    for( unsigned u=1; u<3; u++)
    {
        dg::blas1::transform( q[u], m_temp0, dg::PLUS<double>(-m_p.tbc));
        dg::blas2::symv( m_dx_N, m_temp0, m_dx[u]);
        dg::blas2::symv( m_dy_N, m_temp0, m_dy[u]);
    }
    for( unsigned u=3; u<6; u++)
    {
        dg::blas2::symv( m_dx_U, q[u], m_dx[u]);
        dg::blas2::symv( m_dy_U, q[u], m_dy[u]);
    }
    for( unsigned u=0; u<3; u++)
    {
        dg::blas2::symv( m_dx_P, psi[u], m_dx[6+u]);
        dg::blas2::symv( m_dy_P, psi[u], m_dy[6+u]);
    }
}

template<class Grid, class IMatrix, class Matrix, class Container>
void PerpDynamics<Grid, IMatrix, Matrix, Container>::add_perp_density_dynamics(
    unsigned s,
    const Container& apar,
    const std::array<Container,4>& psi,
    const std::array<std::vector<Container>,6> >& y,
    const std::array<Container,6>& q,
    std::array<std::vector<Container>,6>& yp
)
{
    /////////////////////////////////////////////////////////////////
    double mu = m_p.mu[s], z = m_p.z[s], beta = m_p.beta;
    dg::blas1::subroutine( [mu, z, beta] DG_DEVICE (
            double N,
            double dxFN, double dyFN,
            double dxBN, double dyBN,
            double dxFPperp, double dyFPperp,
            double dxBPperp, double dyBPperp,
            double dxFPpara, double dyFPpara,
            double dxBPpara, double dyBPpara,
            double A, double dxA, double dyA,
            double Tperp, double dxTperp, double dyTperp,
            double Tpara, double dxTpara, double dyTpara,
            double U, double dxU, double dyU,
            double Uperp, double dxUperp, double dyUperp,
            double Upara, double dxUpara, double dyUpara,
            double dxG1, double dyG1,
            double G2, double dxG2, double dyG2,
            double G3, double dxG3, double dyG3,
            double curvNablaX, double curvNablaY,
            double curvKappaX, double curvKappaY,
            double gradLnBX, double gradLnBY,
            double divCurvKappa, double b_2, double divb,
            double& dtN, double& dtPperp, double& dtPpara
        )
    {
        double E0X = dxG1 - G2*(dxTperp/Tperp - gradLnBX);
        double E0Y = dyG1 - G2*(dyTperp/Tperp - gradLnBY);
        double E1X = dxG2 - (G3-G2)*(dxTperp/Tperp - gradLnBX);
        double E1Y = dyG2 - (G3-G2)*(dyTperp/Tperp - gradLnBY);
        double bpX = 0., bpY = 0., divbp = 0.;
        if( beta != 0)
        {
            bpX = A * curvKappaX + (   dyA*b_2);
            bpY = A * curvKappaY + ( - dxA*b_2);
            divbp = A*divCurvKappa - curvNablaX*dxA - curvNablaY*dyA;
        }
        double vX = U * bpX + ( - b_2*E0Y) + Tperp/z *curvNablaX + (Tpara + mu*U*U)/z*curvKappaX;
        double vY = U * bpY + (   b_2*E0X) + Tperp/z *curvNablaY + (Tpara + mu*U*U)/z*curvKappaY;
        dtN += ( vX > 0 ) ? -vX*dxBN : -vX*dxFN;
        dtN += ( vY > 0 ) ? -vY*dyBN : -vY*dyFN;
        // Pperp
        vX = (U+Uperp)*bpX + ( - b_2*(E0Y+E1Y)) + 2*Tperp/z *curvNablaX
            + (Tpara + 2*mu*Uperp*U + mu*U*U)/z*curvKappaX;
        vY = (U+Uperp)*bpY + (   b_2*(E0X+E1X)) + 2*Tperp/z *curvNablaY
            + (Tpara + 2*mu*Uperp*U + mu*U*U)/z*curvKappaY;
        dtPperp += ( vX > 0 ) ? -vX*dxBPperp : -vX*dxFPperp;
        dtPperp += ( vY > 0 ) ? -vY*dyBPperp : -vY*dyFPperp;
        // Ppara
        vX = (U+Upara)*bpX + ( - b_2*E0Y) + Tperp/z *curvNablaX
            + (3*Tpara + 2*mu*Upara*U + mu*U*U)/z*curvKappaX;
        vY = (U+Upara)*bpY + (   b_2*E0X) + Tperp/z *curvNablaY
            + (3*Tpara + 2*mu*Upara*U + mu*U*U)/z*curvKappaY;
        dtPpara += ( vX > 0 ) ? -vX*dxBPpara : -vX*dxFPpara;
        dtPpara += ( vY > 0 ) ? -vY*dyBPpara : -vY*dyFPpara;

        dtN     -= N*( U*divbp + bpX*dxU + bpY*dyU );
        dtPperp -= N*Tperp*( (U+Uperp)*divbp + bpX*(dxU+dxUperp) + bpY*(dyU+dyUperp) );
        dtPpara -= N*Tpara*( (U+Upara)*divbp + bpX*(dxU+dxUpara) + bpY*(dyU+dyUpara) );

        dtN     -= N/z*(Tpara + mu*U*U - Tperp)*divCurvKappa;
        dtPperp -= N*Tperp/z*(Tpara + 2*mu*U*Uperp +  mu*U*U - 2*Tperp)*divCurvKappa;
        dtPpara -= N*Tpara/z*(3*Tpara + 2*mu*U*Upara +  mu*U*U - Tperp)*divCurvKappa;

        dtN     -= N/z*( curvKappaX*(dxTpara + 2*mu*U*dxU) + curvKappaY*(dyTpara + 2*mu*U*dyU));
        dtPperp -= N*Tperp/z*(curvKappaX*(dxTpara + 2*mu*U*dxUperp + 2*mu*(U+Uperp)*dxU)
                              +curvKappaY*(dyTpara + 2*mu*U*dyUperp + 2*mu*(U+Uperp)*dyU));
        dtPpara -= N*Tpara/z*(curvKappaX*(3*dxTpara + 2*mu*U*dxUpara + 2*mu*(U+Upara)*dxU)
                              +curvKappaY*(3*dyTpara + 2*mu*U*dyUpara + 2*mu*(U+Upara)*dyU));

        dtN     -= N/z*(curvNablaX*dxTperp + curvNablaY*dyTperp);
        dtPperp -= N*Tperp/z*(2*curvNablaX*dxTperp + 2*curvNablaY*dyTperp);
        dtPpara -= N*Tpara/z*(curvNablaX*dxTperp + curvNablaY*dyTperp);

        double divuE0 =  (curvNablaX+curvKappaX)*E0X
                        +(curvNablaY+curvKappaY)*E0Y
                        +b_2* ( dxG2*(dyTperp/Tperp - gradLnBY)
                                 -dyG2*(dxTperp/Tperp - gradLnBX));
        double divuE1 =  (curvNablaX+curvKappaX)*E1X
                        +(curvNablaY+curvKappaY)*E1Y
                        +b_2* ( (dxG3-dxG2)*(dyTperp/Tperp - gradLnBY)
                                 -(dyG3-dyG2)*(dxTperp/Tperp - gradLnBX));
        dtN     -= N*divuE0;
        dtPperp -= N*Tperp*(divuE0 + divuE1);
        dtPpara -= N*Tpara*divuE0;


        // F1
        dtPperp -=  N*Tperp*( (divb + divbp)*(Uperp+U)
                   + (Tpara + 2*mu*U*Uperp + mu*U*U)*divCurvKappa/z
                   + curvNablaX*(E0X+E1X) + curvNablaY*(E0Y+E1Y) );
        // F2
        dtPpara += 2*N*Tperp*( (divb + divbp)*Uperp + (Tpara + mu*U*Uperp)*divCurvKappa/z);
        dtPpara -= 2*N*( z*Uperp*(bpX*E1X + bpY*E1Y)
                         + Tpara*(curvKappaX*E0X+curvKappaY*E0Y)
                         + mu*U*Uperp*(curvKappaX*E1X + curvKappaY*E1Y));
        dtPara -=2*N*(bpX*Tpara
                    + mu/z*(Tpara*Upara + 2*U*Tpara)*curvKappaX
                    + mu/z*Tperp*Uperp*curvNablaX
                    -mu*Uperp*b_2*E1Y)*dxU;
        dtPara -=2*N*(bpY*Tpara
                    + mu/z*(Tpara*Upara + 2*U*Tpara)*curvKappaY
                    + mu/z*Tperp*Uperp*curvNablaY
                    +mu*Uperp*b_2*E1X)*dyU;
    },
    y[0][s],
    m_dxF[0], m_dyF[0],
    m_dxB[0], m_dyB[0],
    m_dxF[1], m_dyF[1],
    m_dxB[1], m_dyB[1],
    m_dxF[2], m_dyF[2],
    m_dxB[2], m_dyB[2],
    apar, m_dx[0], m_dy[0],
    q[1], m_dx[1], m_dy[1],
    q[2], m_dx[2], m_dy[2],
    q[3], m_dx[3], m_dy[3],
    q[4], m_dx[4], m_dy[4],
    q[5], m_dx[5], m_dy[5],
    m_dx[6], m_dy[6],
    psi[1], m_dx[7], m_dy[7],
    psi[2], m_dx[8], m_dy[8],
    m_curvNabla[0], m_curvNabla[1],
    m_curvKappa[0], m_curvKappa[1],
    m_gradLnB[0], m_gradLnB[1],
    m_divCurvKappa, m_b_2, m_divb,
    yp[0][s], yp[1][s], yp[2][s]
    );
}

template<class Grid, class IMatrix, class Matrix, class Container>
void PerpDynamics<Grid, IMatrix, Matrix, Container>::add_perp_velocity_dynamics(
    unsigned s,
    const Container& aparST,
    const std::array<Container,4>& psiST,
    const std::array<std::vector<Container>,6> >& y,
    const std::array<Container,6>& qST,
    std::array<std::vector<Container>,6>& yp
)
{
    const Container& to_advect[3] = {qST[3], y[4][s], y[5][s]};
    for( unsigned u=0; u<3; u++)
    {
        dg::blas2::symv( m_dxF_N, to_advect[u], m_dxF[u]);
        dg::blas2::symv( m_dxB_N, to_advect[u], m_dxB[u]);
        dg::blas2::symv( m_dyF_N, to_advect[u], m_dyF[u]);
        dg::blas2::symv( m_dyB_N, to_advect[u], m_dyB[u]);
    }
    for( unsigned u=0; u<3; u++)
    {
        dg::blas1::transform( qST[u], m_temp0, dg::PLUS<double>(u == 0 ? -m_p.nbc[s] : -m_p.tbc));
        dg::blas2::symv( m_dx_N, m_temp0, m_dx[u]);
        dg::blas2::symv( m_dy_N, m_temp0, m_dy[u]);
    }
    dg::blas2::symv( m_dx_A, aparST, m_dx[3]);
    dg::blas2::symv( m_dy_A, aparST, m_dy[3]);
    for( unsigned u=0; u<4; u++)
    {
        dg::blas2::symv( m_dx_P, psiST[u], m_dx[4+u]);
        dg::blas2::symv( m_dy_P, psiST[u], m_dy[4+u]);
    }
    double mu = m_p.mu[s], z = m_p.z[s], beta = m_p.beta;
    dg::blas1::subroutine( [mu, z, beta] DG_DEVICE (
            double U,
            double dxFU, double dyFU,
            double dxBU, double dyBU,
            double dxFQperp, double dyFQperp,
            double dxBQperp, double dyBQperp,
            double dxFQpara, double dyFQpara,
            double dxBQpara, double dyBQpara,
            double N, double dxN, double dyN,
            double Tperp, double dxTperp, double dyTperp,
            double Tpara, double dxTpara, double dyTpara,
            double A, double dxA, double dyA,
            double Uperp, double Upara,
                       double dxG1, double dyG1,
            double G2, double dxG2, double dyG2,
            double G3, double dxG3, double dyG3,
            double G4, double dxG4, double dyG4,
            double curvNablaX, double curvNablaY,
            double curvKappaX, double curvKappaY,
            double gradLnBX, double gradLnBY,
            double divCurvKappa, double b_2, double divb,
            double& dtU, double& dtQperp, double& dtQpara
        )
    {
        double E0X = dxG1 - G2*(dxTperp/Tperp - gradLnBX);
        double E0Y = dyG1 - G2*(dyTperp/Tperp - gradLnBY);
        double E1X = dxG2 - (G3-G2)*(dxTperp/Tperp - gradLnBX);
        double E1Y = dyG2 - (G3-G2)*(dyTperp/Tperp - gradLnBY);
        double E2X = dxG3 - (G4-2*G3)*(dxTperp/Tperp - gradLnBX);
        double E2Y = dyG3 - (G4-2*G3)*(dyTperp/Tperp - gradLnBY);
        double bpX = 0., bpY = 0., divbp = 0.;
        if( beta != 0)
        {
            bpX = A * curvKappaX + (   dyA*b_2);
            bpY = A * curvKappaY + ( - dxA*b_2);
            divbp = A*divCurvKappa - curvNablaX*dxA - curvNablaY*dyA;
        }
        // U
        double vX = U * bpX + ( - b_2*E0Y) + Tperp/z *curvNablaX + (3*Tpara + mu*U*U)/z*curvKappaX;
        double vY = U * bpY + (   b_2*E0X) + Tperp/z *curvNablaY + (3*Tpara + mu*U*U)/z*curvKappaY;
        dtU += ( vX > 0 ) ? -vX*dxBU : -vX*dxFU;
        dtU += ( vY > 0 ) ? -vY*dyBU : -vY*dyFU;
        // Qperp
        vX = U * bpX + ( - b_2*(E0Y+E2Y)) + 3*Tperp/z *curvNablaX + (Tpara + mu*U*U)/z*curvKappaX;
        vY = U * bpY + (   b_2*(E0X+E2X)) + 3*Tperp/z *curvNablaY + (Tpara + mu*U*U)/z*curvKappaY;
        dtQperp += ( vX > 0 ) ? -vX*dxBQperp : -vX*dxFQperp;
        dtQperp += ( vY > 0 ) ? -vY*dyBQperp : -vY*dyFQperp;
        // Qpara
        vX = U * bpX + ( - b_2*E0Y) + Tperp/z *curvNablaX + (7*Tpara + mu*U*U)/z*curvKappaX;
        vY = U * bpY + (   b_2*E0X) + Tperp/z *curvNablaY + (7*Tpara + mu*U*U)/z*curvKappaY;
        dtQpara += ( vX > 0 ) ? -vX*dxBQpara : -vX*dxFQpara;
        dtQpara += ( vY > 0 ) ? -vY*dyBQpara : -vY*dyFQpara;

        dtU -= 2/z*( U*Tpara*divCurvKappa + U * (curvKappaX*dxTpara + curvKappaY*dyTpara)
                        + U*Tpara/N *(curvKappaX*dxN + curvKappaY*dyN))

        double divuE1 =  (curvNablaX+curvKappaX)*E1X
                        +(curvNablaY+curvKappaY)*E1Y
                        +b_2* ( (dxG3-dxG2)*(dyTperp/Tperp - gradLnBY)
                                 -(dyG3-dyG2)*(dxTperp/Tperp - gradLnBX));
        dtU -=    Tpara/mu*(divb+divbp)
                + 1./mu*(bpX*dxTpara + bpY*dyTpara)
                + Tpara/mu/N*(bpX*dxN + bpY*dyN)
                + Tpara*Upara/z*divCurvKappa
                + 1./z/N*( curvKappaX * dxFQpara + curvKappaY * dyFQpara)
                - Tperp/z*Uperp*divCurvKappa
                + 1/z/N*( curvNablaX * dxFQperp + curvNablaY * dyFQperp)
                + Uperp*divuE1
                + 1./N/Tperp*b_2*(E1X*dyFQperp - E1Y*dxFQperp)
                - Uperp/Tperp*b_2*(E1X*dyTperp - E1Y*dxTperp);
        dtU += Tperp/mu*(div+divbp)
                + Tperp/z*(Uperp+U)*divCurvKappa
                - z/mu*(bpX*E0X+bpY*E0Y)
                - U*(    curvKappaX*E0X + curvKappaY*E0Y)
                - Uperp*(curvKappaX*E1X + curvKappaY*E1Y);

        double divuE0 =  (curvNablaX+curvKappaX)*E0X
                        +(curvNablaY+curvKappaY)*E0Y
                        +b_2* ( dxG2*(dyTperp/Tperp - gradLnBY)
                                 -dyG2*(dxTperp/Tperp - gradLnBX));
        double divuE2 =  (curvNablaX+curvKappaX)*E2X
                        +(curvNablaY+curvKappaY)*E2Y
                        +b_2* ( (dxG4-2*dxG3)*(dyTperp/Tperp - gradLnBY)
                                 -(dyG4-2*dyG3)*(dxTperp/Tperp - gradLnBX));
        dtQperp -= N*Tperp*Uperp*(
                        U*divbp + bpX*dxU + bpY*dyU
                        +1/z*(3*Tpara + mu*U*U)*divCurvKappa
                        +1/z*(3*(curvKappaX*dxTpara+curvKappaY*dyTpara)+2*mu*U*(curvKappaX*dxU+curvKappaY*dyU))
                        -3/z*Tperp*divCurvKappa
                        +3/z*(curvNablaX*dxTperp + curvNablaY*dyTperp)
                        +divuE0+divuE2);
        dtQpara -= N*Tpara*Upara*(
                        U*divbp + bpX*dxU + bpY*dyU
                        +1/z*(7*Tpara + mu*U*U)*divCurvKappa
                        +1/z*(7*(curvKappaX*dxTpara+curvKappaY*dyTpara)+2*mu*U*(curvKappaX*dxU+curvKappaY*dyU))
                        -1/z*Tperp*divCurvKappa
                        +1/z*(curvNablaX*dxTperp + curvNablaY*dyTperp)
                        +divuE0);
        // temperature transfer
        vX = N*Tpara/mu*bpX + 1/z*N*Tpara*(Upara+2*mu*U)*curvKappaX + 1/z*N*Tperp*Uperp*curvNablaX
            +N*UPerp*( -b_2*E1Y);
        vY = N*Tpara/mu*bpY + 1/z*N*Tpara*(Upara+2*mu*U)*curvKappaY + 1/z*N*Tperp*Uperp*curvNablaY
            +N*UPerp*(  b_2*E1X);
        dtQperp -= vX*dxTperp + vY*dyTperp;
        dtQpara -= 3*(vX*dxTpara + vY*dyTpara);
        // velocity transfer
        dtQperp -= (N*Tperp*Uperp*bpX + 2*mu/z*U*N*Tperp*Uperp*curvKappaX + 1/z*N*Tperp*Tperp*curvNablaX
            +N*TPerp*( -b_2*E1Y))*dxU;
        dtQperp -= (N*Tperp*Uperp*bpY + 2*mu/z*U*N*Tperp*Uperp*curvKappaY + 1/z*N*Tperp*Tperp*curvNablaY
            +N*TPerp*(  b_2*E1X))*dyU;

        dtQpara -= (N*Tpara*Upara*bpX + 2/z*(N*Tpara*Tpara + mu*U*N*Tpara*Upara)*curvKappaX)*dxU;
        dtQpara -= (N*Tpara*Upara*bpY + 2/z*(N*Tpara*Tpara + mu*U*N*Tpara*Upara)*curvKappaY)*dyU;

        // Force terms
        dtQperp +=   (N*Tperp*(Tperp-Tpara)/mu-N*Tperp*Uperp*U)*(divb+divbp)
                    +1/z * ( 3*N*Tperp*Uperp*(Tperp-Tpara) + N*Tperp*U*(Tperp - 2*Tpara)
                            -Tperp*N*Tpara*Upara - mu*N*Tperp*Uperp*U*U ) * divCurvKappa
                    -(z/mu*N*Tperp*(bpX*E1X+bpY*E1Y)
                        +N*Tperp*Uperp*(curvKappaX*(E0X+E2X)+curvKappaY*(E0Y+E2Y))
                        +N*Tperp*U    *(curvKappaX*E1X + curvKappaY*E1Y)
                        +N*Tperp*Uperp*(curvNablaX*(E0X+E1X+E2X) + curvNablaY*(E0Y+E1Y+E2Y)));
        dtQpara += 3*( (2*Tpara/z*N*Tperp*Uperp + Tperp/z*N*Tpara*Upara) * divCurvKappa
                        -N*Tpara*Upara*(curvKappaX*E0X + curvKappaY*E0Y)
                        -2*N*Tpara*Uperp*(curvKappaX*E1X + curvKappaY*E1Y));

    }
    qST[3],
    m_dxF[0], m_dyF[0],
    m_dxB[0], m_dyB[0],
    m_dxF[1], m_dyF[1],
    m_dxB[1], m_dyB[1],
    m_dxF[2], m_dyF[2],
    m_dxB[2], m_dyB[2],
    qST[0], m_dx[0], m_dy[0],
    qST[1], m_dx[1], m_dy[1],
    qST[2], m_dx[2], m_dy[2],
    aparST, m_dx[3], m_dy[3],
    qST[4], qST[5],
        m_dx[4], m_dy[4],
    psiST[1], m_dx[5], m_dy[5],
    psiST[2], m_dx[6], m_dy[6],
    psiST[3], m_dx[7], m_dy[7],
    m_curvNabla[0], m_curvNabla[1],
    m_curvKappa[0], m_curvKappa[1],
    m_gradLnB[0], m_gradLnB[1],
    m_divCurvKappa, m_b_2, m_divb,
    yp[3][s], yp[4][s], yp[5][s]
    );
}
} //namespace thermal

