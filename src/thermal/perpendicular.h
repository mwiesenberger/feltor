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

    void add_perp_density_dynamics(
        unsigned s,
        const Container& apar,
        const std::array<std::vector<Container>,6> >& y,
        const std::array<Container,6>& q,
        std::array<std::vector<Container>,6>& yp
    );
    private:
    //these should be considered const
    std::array<Container,2> m_curvNabla, m_curvKappa, m_gradLnB;
    Container m_divCurvKappa, m_bphiB, m_divb; //m_bphiB is bhat/ sqrt(g) / B (covariant component)

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
    dg::pushForward(bhat.x(), bhat.y(), bhat.z(), m_temp0, m_temp1, m_bphiB, g);
    // make bhat covariant:
    dg::tensor::inv_multiply3d( g.metric(), m_temp0, m_temp1, m_bphiB,
                                            m_temp0, m_temp1, m_bphiB);
    // Grad Ln B covariant components
    m_gradLnB[0] = dg::pullback( dg::geo::BR(mag), g);
    m_gradLnB[1] = dg::pullback( dg::geo::BZ(mag), g);
    m_temp0 = dg::pullback(dg::geo::InvB(mag), g);
    dg::blas1::pointwiseDot( m_gradLnB[0], m_temp0, m_gradLnB[0]);
    dg::blas1::pointwiseDot( m_gradLnB[1], m_temp0, m_gradLnB[1]);
    dg::blas1::pointwiseDot( m_bphiB, m_temp0, m_bphiB); // bphi / B
    m_temp0 = dg::tensor::volume( metric);
    dg::blas1::pointwiseDivide( m_bphiB, m_temp0, m_bbphiB); //b_2/detg/B

    // allocate resources for internal vectors
    std::fill( m_dxF.begin(), m_dxF.end(), m_temp0);
    m_dxB = m_dyF = m_dyB = m_dxF;
    std::fill( m_dx.begin(), m_dx.end(), m_temp0);
    std::fill( m_dy.begin(), m_dy.end(), m_temp0);
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
    const Container& to_advect[3] = {y[s][0], y[s][1], y[s][2]};
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
    for( unsigned u=6; u<9; u++)
    {
        dg::blas2::symv( m_dx_P, q[u], m_dx[u]);
        dg::blas2::symv( m_dy_P, q[u], m_dy[u]);
    }
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
            double divCurvKappa, double bphiB, double divb,
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
            bpX = A * curvKappaX + (   dyA*bphiB);
            bpY = A * curvKappaY + ( - dxA*bphiB);
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
                        +bphiB* ( dxG2*(dyTperp/Tperp - gradLnBY)
                                 -dyG2*(dxTperp/Tperp - gradLnBX));
        double divuE1 =  (curvNablaX+curvKappaX)*E1X
                        +(curvNablaY+curvKappaY)*E1Y
                        +bphiB* ( (dxG3-dxG2)*(dyTperp/Tperp - gradLnBY)
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
                    -mu*Uperp*bphiB*E1Y)*dxU;
        dtPara -=2*N*(bpY*Tpara
                    + mu/z*(Tpara*Upara + 2*U*Tpara)*curvKappaY
                    + mu/z*Tperp*Uperp*curvNablaY
                    +mu*Uperp*bphiB*E1X)*dyU;
    },
    to_advect[0],
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
    m_divCurvKappa, m_bphiB, m_divb,
    yp[0], yp[1], yp[2]
    );
}

} //namespace thermal

