#pragma once

#include "dg/algorithm.h"
#include "parameters.h"
#include "dg/matrix/matrix.h"
#include "dg/geometries/geometries.h"


namespace thermal
{

// Class to hold solvers for phi and Gamma
// Assume Product Space Grid where phi phi component decouples from perp grids

template< class Geometry, class Matrix, class Container >
struct ThermalSolvers
{
    ThermalSolvers() = default;
    ThermalSolvers( const Geometry&, thermal::Parameters p,
        dg::geo::TokamakMagneticField mag, dg::file::WrappedJsonValue js);
    void compute_phi(
        double time, const std::vector<Container>& density,
        const std::vector<Container>& tperp,
        const std::function<void(Container&)>& multiply_rhs_penalization,
        Container& phi
    );
    void compute_aparST( double t, const std::vector<Container>&,
            std::vector<Container>&, Container&, bool);

    void compute_psi( double time, const Container& phi,
        std::array<Container,3>& psi, unsigned s);

    prviate:
    thermal::Parameters m_p;
    Container m_temp0, m_temp1, m_temp2, m_rhoinv2, m_B, m_vol;
    dg::MultigridCG2d<Geometry, Matrix, Container> m_multigrid;
    std::vector<Container> m_multi_chi;

    std::vector<dg::Elliptic2d< Geometry, Matrix, Container> > m_multi_pol;
    std::vector<dg::Helmholtz2d<Geometry, Matrix, Container> >
        m_multi_invgammaN, m_multi_ampere;

    dg::Elliptic2d<Geometry, Matrix, Container> m_laplaceM;
    dg::mat::ProductMatrixFunction<Container> m_prod;
    cusp::dia_matrix<int, double, cusp::host_memory> m_T;

    dg::MultigridCG2d<Geometry, Matrix, Container> m_multigrid;
    // where do we need _old_apar
    dg::Extrapolation<Container> m_old_phi, m_old_aparST;
    std::vector<dg::Extrapolation<Container>> m_old_gammaN;
};

template<class Geometry, class Matrix, class Container>
Explicit<Geometry, Matrix, Container>::Explicit( const Geometry& g,
    thermal::Parameters p, dg::geo::TokamakMagneticField mag,
    dg::file::WrappedJsonValue js
    ): m_p(p),
    m_multigrid( g, p.stages),
    m_old_phi( 2, dg::evaluate( dg::zero, g)), m_old_aparST( m_old_phi),
    m_old_gammaN( p.num_species - p.num_trivial, m_old_phi),
{
    dg::assign( dg::evaluate( dg::zero, g), m_temp0 );
    m_rhoinv2 = m_temp2, m_temp1 = m_temp0;
    dg::assign(  dg::pullback(dg::geo::Bmodule(mag), g), m_B);
    m_vol = dg::create::volume( g);
    //dg::blas1::pointwiseDot( m_B2, m_B2, m_B2);

    //Set a hard code limit on the maximum number of iteration to avoid
    //endless iteration in case of failure
    m_multigrid.set_max_iter( 1e5);
    /////////////////////////init elliptic and helmholtz operators/////////
    m_multi_chi = m_multigrid.project( m_temp0);
    m_multi_pol.resize(p.stages);
    m_multi_invgammaN.resize(p.stages);
    m_multi_ampere.resize(p.stages);
    for( unsigned u=0; u<p.stages; u++)
    {
        m_multi_pol[u].construct( m_multigrid.grid(u),
            p.bcxP, p.bcyP,
            p.pol_dir, p.jfactor);
        m_multi_invgammaN[u] = { -1.,
                {m_multigrid.grid(u), p.bcxN, p.bcyN, p.pol_dir}};
        m_multi_ampere[u] = {  -1.,
                {m_multigrid.grid(u), p.bcxA, p.bcyA, p.pol_dir}};
    }
    m_laplaceM.construct( g, p.bcxP, p.bcyP, p.pol_dir, p.jfactor);
    m_prod.construct( m_temp0, 1e3);
}

template<class Geometry, class Matrix, class Container>
void ThermalSolvers<Geometry, Matrix, Container>::compute_phi(
    double time, const std::vector<Container>& density,
    const std::vector<Container>& pperp, // bar P_perp
    const std::function<void(Container&)>& multiply_rhs_penalization,
    Container& phi
    )
{
    //----------Compute and set chi----------------------------//
    dg::blas1::copy( 0., m_temp0);
    for( unsigned s = 0; s<m_p.num_species; s++)
    {
        if( !m_p.neglect_mass[s] )
        {
            dg::blas1::pointwiseDivide( m_p.mu[s], density[s], m_B, 0., m_temp1);
            dg::blas1::pointwiseDivide( 1., m_temp1, m_B, 1., m_temp0);
        }
    }
    m_multigrid.project( m_temp0, m_multi_chi);
    for( unsigned u=0; u<m_p.stages; u++)
        m_multi_pol[u].set_chi( m_multi_chi[u]);

    //----------Compute right hand side------------------------//
    dg::blas1::copy( 0, m_temp0);
    double min = 0.;
    for( unsigned s = 0; s<m_p.num_species; s++)
    {
        if( m_p.neglect_mass[s] )
        {
            dg::blas1::transform( density[s], m_temp1, dg::PLUS<double>(-m_p.nbc[s]));
            dg::blas1::axpby( m_p.z[s], m_temp1, 1., m_temp0);
            continue;
        }
        // compute 2/rho_s^2
        dg::blas1::pointwiseDivide( pperp[s], density[s], m_temp1); // bar T_perp = bar P_perp / N
        dg::blas1::pointwiseDivide( 2.*m_p.z[s]*m_p.z[s]/m_p.mu[s], m_B, m_temp1, 0., m_rhoinv2);
        min = std::min( dg::blas1::reduce( m_rhoinv2, 1e308, thrust::minimum<double>()), min);

        m_multigrid.project( m_rhoinv2, m_multi_chi);
        for( unsigned u=0; u<m_p.stages; u++)
            m_multi_invgammaN[u].set_chi( m_multi_chi[u]);

        //compute Gamma^dagger N_s
        dg::blas1::transform( density[s], m_temp1, dg::PLUS<double>(-m_p.nbc[s]));
        m_old_gammaN[s].extrapolate( time, m_temp2);
        m_multigrid.set_benchmark( true, "Gamma N"+m_p.name[s]+"     ");
        std::vector<unsigned> numberG = m_multigrid.solve(
            m_multi_invgammaN, m_temp2, m_temp1, m_p.eps_gamma);
        m_old_gammaN[s].update( time, m_temp2);

        // gamma *2 / rho^2
        dg::blas1::pointwiseDot( m_temp2, m_rhoinv2, m_temp2);

        dg::blas1::axpby( m_p.z[s], m_temp2, 1., m_temp0);
    }
    // Add penalization method
    multiply_rhs_penalization( m_temp0);
    //----------Invert polarisation----------------------------//
    m_old_phi.extrapolate( time, phi);
    m_multigrid.set_benchmark( true, "Polarisation");
    std::vector<unsigned> number = m_multigrid.solve(
        m_multi_pol, phi, m_temp0, m_p.eps_pol);
    m_old_phi.update( time, phi);

    // compute Lanczos tridiagonalisation of Phi
    dg::Timer t;
    t.tic();
    dg::mat::GyrolagK<double> func(0, -1.);
    auto unary_func = dg::mat::make_FuncEigen_Te1( [&](double x) {return func( x, 1./min);});
    m_T = krylovproduct.lanczos().tridiag( unary_func, m_laplaceM, phi, m_vol, m_p.eps_pol, 1.,
                "universal", 1.0, 1);
    t.toc();
#ifdef MPI_VERSION
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
    DG_RANK0 std::cout << "# Lanczos tridiag "<<m_T.iter<<" iterations took "<<t.diff()<<"s\n";



}

template<class Geometry, class IMatrix, class Matrix, class Container>
void ThermalSolvers<Geometry, IMatrix, Matrix, Container>::compute_aparST(
    double time, const std::vector<Container>& densityST,
    const std::vector<Container>& wST, Container& aparST,
    // update determines if old solution is updated (or start iteration from 0)
    bool update)
{
    // beta is nonzero when this function is called
    //----------Compute and set chi----------------------------//
    dg::blas1::copy( 0, m_temp0);
    for( unsigned s=0; s<m_p.num_species; s++)
    {
        dg::blas1::axpby(  m_p.beta*m_p.z[s]*m_p.z[s]/m_p.mu[s],
            m_densityST[s], 1., m_temp0);
    }
    m_multigrid.project( m_temp0, m_multi_chi);
    for( unsigned u=0; u<m_p.stages; u++)
        m_multi_ampere[u].set_chi( m_multi_chi[u]);

    //----------Compute right hand side------------------------//
    dg::blas1::copy( 0, m_temp0);
    for( unsigned s=0; s<m_p.num_species; s++)
    {
        dg::blas1::pointwiseDot(  m_p.beta*m_p.z[s], densityST[s], wST[s],
                                  1., m_temp0);
    }
    //----------Invert Induction Eq----------------------------//
    if( update)
        m_old_aparST.extrapolate( time, aparST);
    m_multigrid.set_benchmark( true, "Apar        ");
    std::vector<unsigned> number = m_multigrid.solve(
        m_multi_ampere, aparST, m_temp0, m_p.eps_ampere);
    if( update)
        m_old_aparST.update( time, aparST);
    if(  number[0] == m_multigrid.max_iter())
        throw dg::Fail( m_p.eps_ampere);
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void ThermalSolvers<Geometry, IMatrix, Matrix, Container>::compute_psi(
    double time, const std::vector<Container>& pperp,
    const Container& phi, std::array<Container,4>& psi, unsigned s
    )
{
    // it's easy to get confused about minusses here
    if( m_p.neglect_mass[s] )
    {
        dg::blas1::copy( phi, psi[0]);
        for( unsigned u=1; u<4; u++)
            dg::blas1::copy( 0, psi[u]);
        return;
    }
    // compute rho_s^2/2
    dg::blas1::pointwiseDivide( pperp[s], density[s], m_temp1); // bar T_perp = bar P_perp / N
    dg::blas1::pointwiseDivide( m_p.mu[s]/2./m_p.z[s]/m_p.z[s], m_temp1, m_B, 0., m_rhoinv2);

    //-----------Solve for Gamma1 Phi---------------------------//
    for( unsigned u=0; u<4; u++)
    {
        dg::mat::GyrolagK<double> func(u, -1.);
        m_prod.compute_vlcl( func, m_rhoinv2, m_laplaceM, m_T, psi[u], phi, m_prod.lanczos().get_bnorm());
    }
    dg::blas1::axpbypgz( -6., m_psi[1], 12., m_psi[2], -6., m_psi[3]);
    dg::blas1::axpby(    -2., m_psi[1],  2., m_psi[2]);
    dg::blas1::scal( -1., m_psi[1]);
    dg::blas1::pointwiseDivide( 1., m_B, m_temp0);
    m_laplaceM.variation( -m_p.mu[s]/2./m_p.z[s], m_temp0, phi, 1., m_psi[0]); // psi_2
}

}//namespace thermal
