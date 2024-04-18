#pragma once

#include "dg/algorithm.h"
#include "parameters.h"
#include "dg/geometries/geometries.h"

#define FELTORPARALLEL 1
#define FELTORPERP 1

namespace thermal
{

struct BPerp{
    //b_perp
    DG_DEVICE void operator()(double A,
        double d0A, double d1A,
        double& bp0, double& bp1, //bperp
        double b_2,
        double curvKappa0,  double curvKappa1
        ){
        bp0 = (   b_2*d1A + A*curvKappa0);
        bp1 = ( - b_2*d0A + A*curvKappa1);
    }
};

template< class Geometry, class IMatrix, class Matrix, class Container >
struct Explicit
{
    // y[0] == density
    // y[1] == velocity
    // y[2] == pperp
    // y[3] == ppara
    using vector = std::array<std::vector<Container>,4>;
    using container = Container;
    Explicit( const Geometry& g, thermal::Parameters p,
        dg::geo::TokamakMagneticField mag, dg::file::WrappedJsonValue js );

    //Given N_i initialize n_e such that phi=0
    void initializene( const Container& ni, Container& ne, std::string initphi);
    //Given n_e initialize N_i such that phi=0
    void initializeni( const Container& ne, Container& ni, std::string initphi);

    void operator()( double t,
        const std::array<std::array<Container,2>,2>& y,
        std::array<std::array<Container,2>,2>& yp);
    void implicit( double t,
        const std::array<std::array<Container,2>,2>& y,
        std::array<std::array<Container,2>,2>& yp);
    void add_implicit_density( double t,
        const std::array<Container,2>& density,
        double beta,
        std::array<Container,2>& yp);
    template<size_t N>
    void add_implicit_velocityST( double t,
        const std::array<Container,2>& densityST,
        const std::array<Container,2>& velocityST,
        double beta,
        std::array<Container,N>& yp);
    /// ///////////////////RESTART    MEMBERS //////////////////////
    const Container& restart_density(int i)const{
        return m_density[i];
    }
    const Container& restart_velocity(int i)const{
        return m_velocityST[i];
    }
    const Container& restart_aparallel() const {
        return m_aparST;
    }
    /// ///////////////////DIAGNOSTIC MEMBERS //////////////////////
    const Geometry& grid() const {
        return m_multigrid.grid(0);
    }
    //potential[0]: electron potential, potential[1]: ion potential
    const Container& uE2() const {
        return m_UE2;
    }
    const Container& density(int i)const{
        return m_density[i];
    }
    const Container&  gammaNi() const{
        if( m_p.tau[1] == 0)
            return m_density[1];
        return m_old_gammaN.head();
    }
    const Container&  gammaPhi() const{
        if( m_p.tau[1] == 0)
            return m_potential[0];
        return m_old_psi.head();
    }


    const Container& density_source(int i)const{
        return m_s[0][i];
    }
    const Container& velocity(int i)const{
        return m_velocity[i];
    }
    const Container& velocity_source(int i){
        update_diag();
        return m_s[1][i];
    }
    const Container& potential(int i) const {
        return m_potential[i];
    }
    const Container& aparallel() const {
        return m_apar;
    }
    const std::array<Container, 3> & gradN (int i) {
        update_diag();
        // m_dFN is updated to the diff_dir direction derivative
        return m_dFN[i];
    }
    const std::array<Container, 3> & gradU (int i) {
        update_diag();
        // m_dFU is updated to the diff_dir direction derivative
        return m_dFU[i];
    }
    const std::array<Container, 3> & gradP (int i) {
        update_diag();
        return m_dP[i];
    }
    const std::array<Container, 3> & gradA () {
        update_diag();
        return m_dA;
    }
    const Container& divNUb( int i) const{
        return m_divNUb[i];
    }
    const Container& dsN (int i) {
        update_diag();
        return m_dsN[i];
    }
    const Container& dsU (int i) {
        update_diag();
        return m_dsU[i];
    }
    const Container& dsP(int i) {
        update_diag();
        return m_dsP[i];
    }
    const Container& dssU(int i){
        update_diag();
        return m_dssU[i];
    }
    const Container& lapParU( unsigned i) {
        update_diag();
        return m_lapParU[i];
    }
    const Container& lapParN( unsigned i) {
        update_diag();
        return m_lapParN[i];
    }
    void compute_gradSN( int i, std::array<Container,3>& gradS) const{
        // MW: don't like this function, if we need more gradients we might
        // want a more flexible solution
        // grad S_ne and grad S_ni
        dg::blas2::symv( m_dxF_N, m_s[0][i], gradS[0]);
        dg::blas2::symv( m_dyF_N, m_s[0][i], gradS[1]);
    }
    void compute_dot_aparallel( Container& tmp) const {
        m_old_apar.derive( tmp);
    }
    const dg::SparseTensor<Container>& projection() const{
        return m_hh;
    }
    const std::array<Container, 2> & curvNabla () const {
        return m_curvNabla;
    }
    const std::array<Container, 2> & curvKappa () const {
        return m_curvKappa;
    }
    const Container& divCurvKappa() const {
        return m_divCurvKappa;
    }
    const Container& bphi( ) const { return m_bphi; }
    const Container& binv( ) const { return m_binv; }
    const Container& divb( ) const { return m_divb; }
    //volume with dG weights
    const Container& vol3d() const { return m_lapperpN.weights();}
    const Container& weights() const { return m_lapperpN.weights();}
    //bhat / sqrt{g} / B (covariant components)
    const Container & bhatgB () const {
        return m_b;
    }
    void compute_lapMperpN (double alpha, const Container& density, Container& temp0, double beta, Container& result)
    {
        // positive Laplacian
        dg::blas1::transform( density, temp0, dg::PLUS<double>(-m_p.nbc));
        dg::blas2::symv( alpha, m_lapperpN, temp0, beta, result);
    }
    void compute_lapMperpU (int i, Container& result)
    {
        dg::blas2::symv( m_lapperpU, m_velocity[i], result);
    }
    void compute_lapMperpP (int i, Container& result)
    {
        m_lapperpP.set_chi( 1.);
        dg::blas2::gemv( m_lapperpP, m_potential[i], result);
    }
    void compute_lapMperpA ( Container& result)
    {
        // only if lapperpU has same direction as lapperpP
        dg::blas2::gemv( m_lapperpU, m_apar, result);
    }
    void compute_bperp( std::array<Container,3>& bperp)
    {
        update_diag();
        dg::blas1::subroutine( BPerp(), m_apar,
            m_dA[0], m_dA[1],
            bperp[0], bperp[1],// bperp on output
            m_b[2],
            m_curvKappa[0], m_curvKappa[1]
        );
    }
    const Container& get_source() const{
        return m_source;
    }
    const Container& get_source_prof() const{
        return m_profne;
    }
    const Container& get_wall() const{
        return m_wall;
    }
    const Container& get_sheath() const{
        return m_sheath;
    }
    const Container& get_sheath_coordinate() const{
        return m_sheath_coordinate;
    }
    void compute_perp_diffusiveN( double alpha, const Container& density,
            Container& temp0, Container& temp1, double beta, Container& result )
    {
        // density = full N
        // result = alpha Lambda_N + beta result
        if( m_p.nu_perp_n > 0)
        {
            dg::blas1::transform( density, temp0, dg::PLUS<double>(-m_p.nbc));
            for( unsigned s=0; s<m_p.diff_order; s++)
            {
                using std::swap;
                swap( temp0, temp1);
                dg::blas2::symv( 1., m_lapperpN, temp1, 0., temp0);
            }
            dg::blas1::axpby( -alpha*m_p.nu_perp_n, temp0, beta, result);
        }
        else
            dg::blas1::scal( result, beta);
    }
    void compute_perp_diffusiveU( double alpha, const Container& velocity,
            const Container& density,
            Container& temp0, Container& temp1, Container& temp2, Container& temp3, double beta, Container& result)
    {
        // density = full N
        // result = alpha Lambda_U + beta result
        if( m_p.nu_perp_u > 0)
        {
            dg::blas1::copy( velocity, temp0);
            for( unsigned s=0; s<m_p.diff_order; s++)
            {
                using std::swap;
                swap( temp0, temp1);
                dg::blas2::symv( 1., m_lapperpU, temp1, 0., temp0);
            }
            if( !m_modify_diff)
                dg::blas1::pointwiseDivide( -alpha*m_p.nu_perp_u, temp0, density, beta, result);
            else
                dg::blas1::axpby( -alpha*m_p.nu_perp_u, temp0, beta, result);
        }
        else
            dg::blas1::scal( result, beta);
        double nu = m_p.nu_perp_n;
        if( m_modify_diff)
            nu += m_p.nu_perp_u;
        if( nu > 0 )
        {

            dg::blas1::transform( density, temp0, dg::PLUS<double>(-m_p.nbc));
            for( unsigned s=0; s<m_p.diff_order-1; s++)
            {
                using std::swap;
                swap( temp0, temp1);
                dg::blas2::symv( 1., m_lapperpN, temp1, 0., temp0);
            }

            // - v_x dx U
            if( m_p.diff_dir == dg::centered)
                dg::blas2::symv( m_dxC, temp0, temp1);
            else if( m_p.diff_dir == dg::forward)
                dg::blas2::symv( m_dxF_N, temp0, temp1);
            else
                dg::blas2::symv( m_dxB_N, temp0, temp1);
            dg::blas1::pointwiseDivide( -nu, temp1, density, 0., temp1);
            dg::blas2::symv( m_dxB_U, velocity, temp2);
            dg::blas2::symv( m_dxF_U, velocity, temp3);
            dg::blas1::evaluate( result, dg::minus_equals(), dg::UpwindProduct(),
                    temp1, temp2, temp3);
            // - v_y dy U
            if( m_p.diff_dir == dg::centered)
                dg::blas2::symv( m_dyC, temp0, temp1);
            else if( m_p.diff_dir == dg::forward)
                dg::blas2::symv( m_dyF_N, temp0, temp1);
            else
                dg::blas2::symv( m_dyB_N, temp0, temp1);
            dg::blas1::pointwiseDivide( -nu, temp1, density, 0., temp1);
            dg::blas2::symv( m_dyB_U, velocity, temp2);
            dg::blas2::symv( m_dyF_U, velocity, temp3);
            dg::blas1::evaluate( result, dg::minus_equals(), dg::UpwindProduct(),
                    temp1, temp2, temp3);
        }
    }

    bool modify_diff() const {return m_modify_diff;}
    void compute_parallel_diffusiveN( int i, Container& result)
    {
        dg::blas1::axpby( m_p.nu_parallel_n, lapParN(i), 0., result);
    }
    void compute_parallel_diffusiveU( int i, Container& result)
    {
        double nu = m_p.nu_parallel_n;
        if( m_modify_diff)
            nu += m_p.nu_parallel_u[i];
        if( nu > 0)
        {
            dg::blas1::pointwiseDot( dsN(i), dsU(i), result);
            dg::blas1::pointwiseDivide( nu, result, density(1), 0., result);
        }
        else
            dg::blas1::copy( 0, result);
        if( m_p.nu_parallel_u[i] > 0)
        {
            if( !m_modify_diff)
                dg::blas1::pointwiseDivide( m_p.nu_parallel_u[i], lapParU(i), density(i), 1., result);
            else
                dg::blas1::axpby( m_p.nu_parallel_u[i], lapParU(i), 1., result);
        }
    }


    // Compute divergence using centered derivatives
    // note that no matter how divergence is computed you always loose one order
    // unless the polarisation term or the Laplacian of N,U is computed
    // then the correct direction must be chosen
    // prefactor cannot alias result!!
    // Div ( f v)
    template<class Container2>
    void centered_div( const Container2& prefactor,
            const std::array<Container, 3>& contra_vec,
            Container& temp0, Container& result)
    {
        dg::blas1::pointwiseDot( 1., prefactor, m_detg, contra_vec[0], 0., temp0);
        dg::blas2::symv( m_dxC, temp0, result);
        dg::blas1::pointwiseDot( 1., prefactor, m_detg, contra_vec[1], 0., temp0);
        dg::blas2::symv( 1., m_dyC, temp0, 1., result);
        dg::blas1::pointwiseDivide( 1., result, m_detg, 0., result);
    }
    void centered_v_dot_nabla( const std::array<Container, 3>& contra_vec,
            const Container& f, Container& temp1, Container& result)
    {
        dg::blas2::symv( m_dxC, f, temp1);
        dg::blas1::pointwiseDot( contra_vec[0], temp1, result);
        dg::blas2::symv( m_dyC, f, temp1);
        dg::blas1::pointwiseDot( 1., contra_vec[1], temp1, 1., result);
    }
    void compute_pol( double alpha, const Container& density, Container& temp, double beta, Container& result)
    {
        // polarisation term
        dg::blas1::pointwiseDot( m_p.mu[1], density, m_binv, m_binv, 0., temp);
        m_multi_pol[0].set_chi( temp);
        dg::blas2::symv( -alpha, m_multi_pol[0], m_potential[0], beta, result);
    }
    void compute_source_pol( double alpha, const Container& density, Container& temp, double beta, Container& result)
    {
        // we don't want jumps in phi in here so we use lapperpP
        dg::blas1::pointwiseDot( m_p.mu[1], density, m_binv, m_binv, 0., temp);
        m_lapperpP.set_chi( temp);
        dg::blas2::symv( -alpha, m_lapperpP, m_potential[0], beta, result);
    }
    unsigned called() const { return m_called;}

    /// //////////////////////DIAGNOSTICS END////////////////////////////////
    void update_diag(){
        // assume m_density, m_potential, m_velocity, m_velocityST, m_apar
        // compute dsN, dsU, dsP, lapParU, dssU and perp derivatives
        if( !m_upToDate)
        {
            // update m_dN, m_dU, m_dP, m_dA
            update_perp_derivatives( m_density, m_velocity, m_potential, m_apar);
            for( unsigned i=0; i<2; i++)
            {
                // density m_dsN, m_lapParN
                m_fa( dg::geo::einsMinus, m_density[i], m_minus);
                m_fa( dg::geo::zeroForw,  m_density[i], m_zero);
                m_fa( dg::geo::einsPlus,  m_density[i], m_plus);
                update_parallel_bc_2nd( m_fa, m_minus, m_zero, m_plus,
                        m_p.bcxN, m_p.bcxN == dg::DIR ? m_p.nbc : 0.);
                dg::geo::ds_centered( m_fa, 1., m_minus, m_plus, 0., m_dsN[i]);
                dg::geo::dssd_centered( m_fa, 1.,
                        m_minus, m_zero, m_plus, 0., m_lapParN[i]);
                // potential m_dsP
                m_fa( dg::geo::einsMinus, m_potential[i], m_minus);
                m_fa( dg::geo::einsPlus,  m_potential[i], m_plus);
                update_parallel_bc_2nd( m_fa, m_minus, m_potential[i], m_plus,
                        m_p.bcxP, 0.);
                dg::geo::ds_centered( m_fa, 1., m_minus, m_plus, 0., m_dsP[i]);
                // velocity m_dssU, m_lapParU m_dsU
                m_fa( dg::geo::einsMinus, m_velocity[i], m_minus);
                m_fa( dg::geo::zeroForw,  m_velocity[i], m_zero);
                m_fa( dg::geo::einsPlus,  m_velocity[i], m_plus);
                update_parallel_bc_2nd( m_fa, m_minus, m_zero, m_plus,
                        m_p.bcxU, 0.);
                dg::geo::dssd_centered( m_fa, 1.,
                        m_minus, m_zero, m_plus, 0., m_lapParU[i]);
                dg::geo::dss_centered( m_fa, 1., m_minus,
                    m_zero, m_plus, 0., m_dssU[i]);
                dg::geo::ds_centered( m_fa, 1., m_minus, m_plus, 0.,
                        m_dsU[i]);
                // velocity source
                dg::blas1::evaluate( m_s[1][i], dg::equals(), []DG_DEVICE(
                            double sn, double u, double n){ return -u*sn/n;},
                        m_s[0][i], m_velocity[i], m_density[i]);
            }
            for( unsigned i=0; i<2; i++)
            for( unsigned j=0; j<3; j++)
            {
                // update m_dFU, m_dFN to the diff_dir direction derivative
                if( m_p.diff_dir == dg::forward)
                {
                    ;
                }
                else if( m_p.diff_dir == dg::backward)
                {
                    dg::blas1::copy( m_dBN[i][j], m_dFN[i][j]);
                    dg::blas1::copy( m_dBU[i][j], m_dFU[i][j]);
                }
                else
                {
                    dg::blas1::axpby( 1./2.,m_dBN[i][j], 1./2., m_dFN[i][j]);
                    dg::blas1::axpby( 1./2.,m_dBU[i][j], 1./2., m_dFU[i][j]);
                }
            }
            m_upToDate = true;
        }

    }
    void update_parallel_bc_1st( Container& minusST, Container& plusST,
            dg::bc bcx, double value)
    {
        if( m_p.fci_bc == "along_field")
            dg::geo::assign_bc_along_field_1st( m_faST, minusST, plusST,
                    minusST, plusST, bcx, {value,value});
        else
        {
            if( bcx == dg::DIR)
            {
                dg::blas1::plus( minusST, -value);
                dg::geo::swap_bc_perp( m_fa, minusST, plusST,
                        minusST, plusST);
                dg::blas1::plus( minusST, +value);
            }
        }
    }
    void update_parallel_bc_2nd( const dg::geo::Fieldaligned<Geometry, IMatrix,
            Container>& fa, Container& minus, const Container& value0,
            Container& plus, dg::bc bcx, double value)
    {
        if( m_p.fci_bc == "along_field")
        {
            dg::geo::assign_bc_along_field_2nd( fa, minus, value0,
                    plus, minus, plus, bcx, {value,value});
        }
        else
        {
            if( bcx == dg::DIR)
            {
                dg::blas1::plus( minus, -value);
                dg::geo::swap_bc_perp( fa, minus, plus,
                        minus, plus);
                dg::blas1::plus( minus, +value);
            }
        }
    }

    //source strength, profile - 1
    void set_source( bool fixed_profile, Container profile, double source_rate, Container source, double minne, double minrate, double minalpha)
    {
        m_fixed_profile = fixed_profile;
        m_profne = profile;
        m_source_rate = source_rate;
        m_source = source;
        m_minne = minne;
        m_minrate = minrate;
        m_minalpha = minalpha;
    }
    void set_wall(double wall_rate, const Container& wall, double nwall, double uwall)
    {
        m_wall_rate = wall_rate;
        dg::blas1::copy( wall, m_wall);
        m_nwall = nwall;
        m_uwall = uwall;
    }
    void set_sheath(double sheath_rate, const Container& sheath,
            const Container& sheath_coordinate)
    {
        m_sheath_rate = sheath_rate;
        dg::blas1::copy( sheath, m_sheath);
        dg::blas1::copy( sheath_coordinate, m_sheath_coordinate);
    }
    void update_staggered_density_and_phi( double t,
        const std::array<Container,2>& density,
        const std::array<Container,2>& potential);
    void update_staggered_density( double t,
        const std::array<Container,2>& density);
    void update_velocity( double t,
        const std::array<Container,2>& velocityST,
        const Container& aparST);

    void update_perp_derivatives(
        const std::array<Container,2>& density,
        const std::array<Container,2>& velocity,
        const std::array<Container,2>& potential,
        const Container& apar);
    void compute_perp_density( double t,
        const std::array<Container,2>& density,
        const std::array<Container,2>& velocity,
        const std::array<Container,2>& potential,
        const Container& apar,
        std::array<Container,2>& densityDOT);
    void compute_perp_velocity( double t,
        const std::array<Container,2>& density,
        const std::array<Container,2>& velocity,
        const std::array<Container,2>& potential,
        const Container& apar,
        std::array<Container,2>& velocityDOT);
    void compute_parallel_flux(
             const Container& velocityKM,
             const Container& velocityKP,
             const Container& densityM,
             const Container& density,
             const Container& densityP,
             Container& fluxM,
             Container& fluxP,
             std::string slope_limiter);
    void compute_parallel_advection(
             const Container& velocityKM,
             const Container& velocityKP,
             const Container& densityM,
             const Container& density,
             const Container& densityP,
             Container& fluxM,
             Container& fluxP,
             std::string slope_limiter);
    void compute_parallel_flux(
        const Container& velocity,
        const Container& minusST,
        const Container& plusST,
        Container& flux,
        std::string slope_limiter);
    void compute_parallel_advection(
        const Container& velocity,
        const Container& minusST,
        const Container& plusST,
        Container& flux,
        std::string slope_limiter);
    void compute_parallel(          std::array<std::array<Container,2>,2>& yp);
    void multiply_rhs_penalization(      Container& yp);
    void add_wall_and_sheath_terms( std::array<std::array<Container,2>,2>& yp);
    void add_source_terms(          std::array<std::array<Container,2>,2>& yp);
    const dg::geo::Fieldaligned<Geometry, IMatrix, Container>& fieldaligned() const
    {
        return m_fa;
    }
  private:
    void construct_mag( const Geometry&, thermal::Parameters,
        dg::geo::TokamakMagneticField);
    void construct_bhat( const Geometry&, thermal::Parameters,
        dg::geo::TokamakMagneticField);

    //these should be considered const // m_curv is full curvature
    std::array<Container,2> m_curvNabla, m_curvKappa;
    Container m_divCurvKappa, m_b; //m_b is bhat/ sqrt(g) / B (covariant component)
    Container m_bphi, m_binv, m_divb, m_detg; // m_bphi is covariant bhat component
    Container m_source, m_profne, m_sheath_coordinate;
    Container m_wall, m_sheath;

    // Only set once every call to operator()
    std::array<Container,2> m_velocity, m_velocityST;
    std::array<Container,2> m_potential, m_potentialST;
    Container m_phi, m_apar, m_aparST;

    Container m_UE2;
    Container m_divNUb;
    Container m_plusN, m_zeroN, m_minusN, m_plusU, m_zeroU, m_minusU;
    std::vector<Container> m_plusSTN, m_minusSTN, m_densityST;
    Container m_plusSTU, m_minusSTU;

    // overwritten by diag_update and/or set once by operator()
    std::array<Container,2> m_dA;
    std::array<Container,2> m_dP, m_dFN, m_dBN, m_dFU, m_dBU;
    Container m_dsN, m_dsP, m_dsU;
    std::array<std::vector<Container>,4> m_s; // source

    // Set by diag_update
    Container m_dssU, m_lapParU, m_lapParN;

    // Helper variables can be overwritten any time (except by compute_parallel)!!
    Container m_temp0, m_temp1;
    Container m_minus, m_zero, m_plus;
    // Helper variables for compute_parallel_flux
    Container m_vbm, m_vbp, m_dN, m_dNMM, m_dNM, m_dNZ, m_dNP, m_dNPP;

    //matrices and solvers
    Matrix m_dxF_N, m_dxB_N, m_dxF_U, m_dxB_U, m_dx_P, m_dx_A;
    Matrix m_dyF_N, m_dyB_N, m_dyF_U, m_dyB_U, m_dy_P, m_dy_A;
    Matrix m_dxC, m_dyC;
    dg::geo::Fieldaligned<Geometry, IMatrix, Container> m_fa, m_faST;
    dg::Elliptic2d< Geometry, Matrix, Container> m_lapperpN, m_lapperpU, m_lapperpP;
    thermal::ThermalSolvers<Geometry, Matrix, Container> m_solvers;

    const thermal::Parameters m_p;
    const dg::file::WrappedJsonValue m_js;
    double m_source_rate = 0., m_sheath_rate = 0., m_wall_rate = 0.;
    double m_minne = 0., m_minrate  = 0., m_minalpha = 0.;
    double m_nwall = 0., m_uwall = 0.;
    bool m_fixed_profile = true, m_reversed_field = false;
    std::vector<bool> m_upToDate;
    bool m_modify_diff = false;
    unsigned m_called = 0;

};

template<class Grid, class IMatrix, class Matrix, class Container>
void Explicit<Grid, IMatrix, Matrix, Container>::construct_mag(
    const Grid& g, thermal::Parameters p, dg::geo::TokamakMagneticField mag
    )
{
    //due to the various approximations bhat and mag not always correspond
    dg::geo::CylindricalVectorLvl0 curvNabla, curvKappa;
    m_reversed_field = false;
    if( mag.ipol()( g.x0(), g.y0()) < 0)
        m_reversed_field = true;
    if( p.curvmode == "true" )
        throw std::runtime_error( "curvmode : true is not possible in thermal code!");
    else if( p.curvmode == "low beta")
    {
        if( m_reversed_field)
            curvNabla = curvKappa = dg::geo::createCurvatureNablaB(mag, -1);
        else
            curvNabla = curvKappa = dg::geo::createCurvatureNablaB(mag, +1);
        dg::assign( dg::evaluate(dg::zero, g), m_divCurvKappa);
    }
    else if( p.curvmode == "toroidal")
    {
        if( m_reversed_field)
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
    dg::assign(  dg::pullback(dg::geo::InvB(mag), g), m_binv);
    dg::assign(  dg::pullback(dg::geo::Divb(mag), g), m_divb);

}
template<class Grid, class IMatrix, class Matrix, class Container>
void Explicit<Grid, IMatrix, Matrix, Container>::construct_bhat(
    const Grid& g, thermal::Parameters p, dg::geo::TokamakMagneticField mag)
{
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

    // in Poisson we take EPhi except for the true curvmode
    bhat = dg::geo::createEPhi(+1);
    if( p.curvmode == "true")
        throw std::runtime_error( "curvmode : true is not possible in thermal code!");
    else if( m_reversed_field)
        bhat = dg::geo::createEPhi(-1);
    dg::pushForward(bhat.x(), bhat.y(), bhat.z(), m_temp0, m_temp1, m_b, g);
    dg::SparseTensor<Container> metric = g.metric();
    // make bhat covariant:
    dg::tensor::inv_multiply3d( metric, m_temp0, m_temp1, m_b,
                                        m_temp0, m_temp1, m_b);
    dg::assign( m_b, m_bphi); //save bphi for momentum conservation
    m_detg = dg::tensor::volume( metric);
    dg::blas1::pointwiseDivide( m_binv, m_detg, m_temp0); //1/B/detg
    dg::blas1::pointwiseDot( m_temp0, m_b, m_b); //b_2/detg/B
    m_lapperpN.construct ( g, p.bcxN, p.bcyN, dg::PER,  p.diff_dir),
    m_lapperpU.construct ( g, p.bcxU, p.bcyU, dg::PER,  p.diff_dir),
    m_lapperpP.construct ( g, p.bcxP, p.bcyP, dg::PER,  p.pol_dir),
    m_lapperpP.set_jfactor(0); //we don't want jump terms in source
}

template<class Grid, class IMatrix, class Matrix, class Container>
Explicit<Grid, IMatrix, Matrix, Container>::Explicit( const Grid& g,
    thermal::Parameters p, dg::geo::TokamakMagneticField mag,
    dg::file::WrappedJsonValue js
    ):
    m_dxF_N( dg::create::dx( g, p.bcxN, dg::forward) ),
    m_dxB_N( dg::create::dx( g, p.bcxN, dg::backward) ),
    m_dxF_U( dg::create::dx( g, p.bcxU, dg::forward) ),
    m_dxB_U( dg::create::dx( g, p.bcxU, dg::backward) ),
    m_dx_P(  dg::create::dx( g, p.bcxP, p.pol_dir) ),
    m_dx_A(  dg::create::dx( g, p.bcxA, p.pol_dir) ),
    m_dyF_N( dg::create::dy( g, p.bcyN, dg::forward) ),
    m_dyB_N( dg::create::dy( g, p.bcyN, dg::backward) ),
    m_dyF_U( dg::create::dy( g, p.bcyU, dg::forward) ),
    m_dyB_U( dg::create::dy( g, p.bcyU, dg::backward) ),
    m_dy_P(  dg::create::dy( g, p.bcyP, p.pol_dir) ),
    m_dy_A(  dg::create::dy( g, p.bcyA, p.pol_dir) ),
    m_dxC(   dg::create::dx( g, dg::NEU, dg::centered) ), // for divergence
    m_dyC(   dg::create::dy( g, dg::NEU, dg::centered) ), // for divergence
    m_old_apar( 2, dg::evaluate( dg::zero, g)),
    m_solvers( g, p, mag, js),
    m_p(p), m_js(js)
{
    //--------------------------init vectors to 0-----------------//
    dg::assign( dg::evaluate( dg::zero, g), m_temp0 );
    m_source = m_sheath_coordinate = m_UE2 = m_temp1 = m_temp0;
    m_apar = m_aparST = m_profne = m_wall = m_sheath = m_temp0;
    m_plus = m_zero = m_minus = m_temp0;
    m_vbm = m_vbp = m_temp0;
    if( m_p.slope_limiter != "none")
        m_dN = m_dNMM = m_dNM = m_dNZ = m_dNP = m_dNPP = m_temp1;

    m_potential[0] = m_potential[1] = m_temp0;
    m_plusSTN = m_minusSTN = m_minusSTU = m_plusSTU = m_potential;
    m_plusN = m_zeroN = m_minusN = m_minusU = m_zeroU = m_plusU = m_potential;
    m_divNUb = m_density = m_densityST = m_velocity = m_potential;
    m_velocityST = m_potentialST = m_potential;
    m_dsN = m_dsP = m_dsU = m_dssU = m_lapParU = m_lapParN = m_potential;

    m_dA[0] = m_dA[1] = m_dA[2] = m_temp0;
    m_dP[0] = m_dP[1] = m_dA;
    m_dFN = m_dBN = m_dFU = m_dBU = m_dP;
    m_s[0] = m_s[1] = m_potential ;

    //--------------------------Construct-------------------------//
    construct_mag( g, p, mag);
    construct_bhat( g, p, mag);
    // An optional hidden parameter
    if( m_js["regularization"].isMember("modify-diff") )
        m_modify_diff = m_js["regularization"]["modify-diff"].asBool();
#ifdef MPI_VERSION
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
    if( m_modify_diff)
        DG_RANK0 std::cout << "# Optional parameter \"modify-diff\" activated\n";
    m_upToDate.resize( m_p.num_species);
    for( unsigned s=0; s<m_p.num_species; s++)
        m_upToDate[s] = false;

}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::initializene(
    const Container& src, Container& target, std::string initphi)
{
    // ne  = Ni
    dg::blas1::copy( src, target);
    if (m_p.tau[1] != 0.) {
        if( initphi == "zero")
        {
            // ne-nbc = Gamma (ni-nbc)
            dg::blas1::transform(src, m_temp0, dg::PLUS<double>(-m_p.nbc));
            dg::blas1::plus(target, -m_p.nbc);
            m_multigrid.set_benchmark( true, "Gamma N     ");
            std::vector<unsigned> number = m_multigrid.solve(
                m_multi_invgammaN, target, m_temp0, m_p.eps_gamma);
            dg::blas1::plus(target, +m_p.nbc);
        }
        else if( initphi == "balance")
        {
            //add FLR correction -0.5*tau*mu*Delta n_e
            dg::blas1::transform(src, m_temp0, dg::PLUS<double>(-m_p.nbc));
            dg::blas2::symv( 0.5*m_p.tau[1]*m_p.mu[1],
                m_lapperpN, m_temp0, 1.0, target);
            //wird stark negativ falls alpha klein!!
        }
        else if( !(initphi == "zero_pol"))
        {
            throw dg::Error(dg::Message(_ping_)<<"Warning! initphi value '"<<initphi<<"' not recognized. I have tau = "<<m_p.tau[1]<<" ! I don't know what to do! I exit!\n");
        }
    }
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::initializeni(
    const Container& src, Container& target, std::string initphi)
{
    //According to Markus we should actually always invert
    //so we should reconsider this function
    // Ni = ne
    dg::blas1::copy( src, target);
    if (m_p.tau[1] != 0.) {
        if( initphi == "zero")
        {
            //add FLR correction -0.5*tau*mu*Delta n_e
            dg::blas1::transform(src, m_temp0, dg::PLUS<double>(-m_p.nbc));
            dg::blas2::symv( 0.5*m_p.tau[1]*m_p.mu[1],
                m_lapperpN, m_temp0, 1.0, target);
            //wird stark negativ falls alpha klein!!
        }
        else if( initphi == "balance")
        {
            //add FLR correction +0.5*tau*mu*Delta n_e
            dg::blas1::transform(src, m_temp0, dg::PLUS<double>(-m_p.nbc));
            dg::blas2::symv( -0.5*m_p.tau[1]*m_p.mu[1],
                m_lapperpN, m_temp0, 1.0, target);
            //wird stark negativ falls alpha klein!!
        }
        else if( !(initphi == "zero_pol"))
        {
            throw dg::Error(dg::Message(_ping_)<<"Warning! initphi value '"<<initphi<<"' not recognized. I have tau = "<<m_p.tau[1]<<" ! I don't know what to do! I exit!\n");
        }
    }
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::update_perp_derivatives(
    const Container& density,
    const Container& velocity,
    const Container& potential,
    const Container& apar)
{
    for( unsigned i=0; i<2; i++)
    {
        ////////////////////perpendicular dynamics////////////////////////
        //First compute forward and backward derivatives for upwind scheme
        dg::blas1::transform( density, m_temp1, dg::PLUS<double>(-m_p.nbc));
        dg::blas2::symv( m_dxF_N, m_temp1, m_dFN[0]);
        dg::blas2::symv( m_dyF_N, m_temp1, m_dFN[1]);
        dg::blas2::symv( m_dxB_N, m_temp1, m_dBN[0]);
        dg::blas2::symv( m_dyB_N, m_temp1, m_dBN[1]);
        dg::blas2::symv( m_dxF_U, velocity, m_dFU[0]);
        dg::blas2::symv( m_dyF_U, velocity, m_dFU[1]);
        dg::blas2::symv( m_dxB_U, velocity, m_dBU[0]);
        dg::blas2::symv( m_dyB_U, velocity, m_dBU[1]);
        dg::blas2::symv( m_dx_P, potential, m_dP[0]);
        dg::blas2::symv( m_dy_P, potential, m_dP[1]);
        dg::blas2::symv( m_dx_A, apar, m_dA[0]);
        dg::blas2::symv( m_dy_A, apar, m_dA[1]);
    }
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::update_staggered_density_and_phi(
    double t,
    const Container& density,
    const Container& potential, unsigned s)
{
    for( unsigned i=0; i<2; i++)
    {
        m_fa( dg::geo::einsMinus, density, m_minusN);
        m_fa( dg::geo::zeroForw,  density, m_zeroN);
        m_fa( dg::geo::einsPlus,  density, m_plusN);
        update_parallel_bc_2nd( m_fa, m_minusN, m_zeroN, m_plusN,
                m_p.bcxN, m_p.bcxN == dg::DIR ? m_p.nbc : 0.);

        m_faST( dg::geo::zeroMinus, potential[i], m_minus);
        m_faST( dg::geo::einsPlus,  potential[i], m_plus);
        update_parallel_bc_1st( m_minus, m_plus,
                m_p.bcxP, 0.);
        dg::geo::ds_centered( m_faST, 1., m_minus, m_plus, 0., m_dsP[i]);
        dg::blas1::axpby( 0.5, m_minus, 0.5, m_plus, m_potentialST[i]);
    }
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::update_staggered_density(
    double t,
    const std::vector<Container>& density)
{
    for( unsigned s=0; s<m_p.num_species; s++)
    {
        m_faST( dg::geo::zeroMinus, density[s], m_minusSTN[s]);
        m_faST( dg::geo::einsPlus,  density[s], m_plusSTN[s]);
        update_parallel_bc_1st( m_minusSTN[s], m_plusSTN[s],
                m_p.bcxN, m_p.bcxN == dg::DIR ? m_p.nbc : 0.);
        dg::blas1::axpby( 0.5, m_minusSTN[s], 0.5, m_plusSTN[s], m_densityST[s]);
    }
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::update_velocity(
    double t,
    const Container& wST,
    const Container& aparST,
    Container& velocityST, unsigned s)
{
    dg::blas1::axpby( 1., wST, -m_p.z[s]/m_p.mu[s], m_aparST, velocityST);
    // Compute dsU and velocity
    m_faST( dg::geo::einsMinus, velocityST, m_minusSTU);
    m_faST( dg::geo::zeroPlus,  velocityST,  m_plusSTU);
    update_parallel_bc_1st( m_minusSTU, m_plusSTU, m_p.bcxU, 0.);
    dg::blas1::axpby( 0.5, m_minusSTU, 0.5, m_plusSTU, m_velocity);

    m_fa( dg::geo::einsMinus, velocityST, m_minusU);
    m_fa( dg::geo::zeroForw,  velocityST, m_zeroU);
    m_fa( dg::geo::einsPlus,  velocityST, m_plusU);
    update_parallel_bc_2nd( m_fa, m_minusU, m_zeroU,
            m_plusU, m_p.bcxU, 0.);
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::compute_perp_density(
    double t,
    const Container& density,
    const Container& velocity,
    const Container& pperp,
    const Container& ppara,
    const std::array<Container,2>& E1,
    const Container& apar,
    double mu, double beta, double z,
    Container& densityDOT)
{
    update_perp_derivatives( density, velocity, potential, apar);
    //y[0] = N, y[1] = W; fields[0] = N, fields[1] = U
    ////////////////////perpendicular dynamics////////////////////////
    dg::direction diff_dir = m_p.diff_dir;
    dg::blas1::subroutine( [mu, z, beta, diff_dir] DG_DEVICE (
            double N, double d0FN, double d1FN,
                      double d0BN, double d1BN,
            double U, double d0FU, double d1FU,
                      double d0BU, double d1BU,
            double Tperp, double d0Tperp, double d1Tperp,
            double Tpara, double d0Tpara, double d1Tpara,
                      double E1R, double E1Z,
            double A, double d0A, double d1A,
            double b_2,
            double curvNabla0, double curvNabla1,
            double curvKappa0, double curvKappa1, double divCurvKappa,
            double& dtN
        )
        {
            dtN = 0;
            // density - upwind scheme
            double v0 = ( - b_2*E1Z) + Tperp/z *curvNabla0 + (Tpara + mu*U*U)/z*curvKappa0;
            double v1 = (   b_2*E1R) + Tperp/z *curvNabla1 + (Tpara + mu*U*U)/z*curvKappa1;
            double bp0 = 0., bp1 = 0., bp2 = 0.;
            if( beta != 0)
            {
                bp0 = A * curvKappa0 + (   d1A*b_2);
                bp1 = A * curvKappa1 + ( - d0A*b_2);

                v0 += U * bp0;
                v1 += U * bp1;
            }
            dtN += ( v0 > 0 ) ? -v0*d0BN : -v0*d0FN;
            dtN += ( v1 > 0 ) ? -v1*d1BN : -v1*d1FN;

            double d0U = 0, d1U = 0, d2U = 0;
            if( diff_dir == dg::forward)
                d0U = d0FU, d1U = d1FU, d2U = d2FU;
            else if( diff_dir == dg::backward)
                d0U = d0BU, d1U = d1BU, d2U = d2BU;
            else
                d0U = (d0FU+d0BU)/2., d1U = (d1FU+d1BU)/2.,
                    d2U = (d2FU+d2BU) / 2.;

            double KappaU = curvKappa0*d0U+curvKappa1*d1U+curvKappa2*d2U;
            double KP =     curvNabla0*d0P+curvNabla1*d1P+curvNabla2*d2P;

            dtN +=  - N * ( KP + mu * U * U * divCurvKappa
                            + 2. * mu * U * KappaU);
            if( beta != 0)
            {
                double divbp = A*divCurvKappa
                                 - curvNabla0*d0A
                                 - curvNabla1*d1A;
                double bpU = bp0*d0U + bp1*d1U;
                dtN +=  -N*( U*divbp + bpU);
            }
            return;
        },
        //species depdendent
        density,   m_dFN[0], m_dFN[1], m_dBN[0], m_dBN[1],
        velocity,  m_dFU[0], m_dFU[1], m_dBU[0], m_dBU[1],
                   m_E1[0], m_E1[1],
        //aparallel
        apar, m_dA[0], m_dA[1],
        //magnetic parameters
        m_b[2], m_curvNabla[0], m_curvNabla[1], m_curvKappa[0], m_curvKappa[1],
        m_divCurvKappa, densityDOT
    );
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::compute_perp_velocity(
    double t,
    const Container& density,
    const Container& velocity,
    const Container& potential,
    const Container& apar,
    double mu, double beta, double z,
    Container& velocityDOT)
{
    update_perp_derivatives( density, velocity, potential, apar);
    //y[0] = N, y[1] = W; fields[0] = N, fields[1] = U
    for( unsigned i=0; i<2; i++)
    {
        ////////////////////perpendicular dynamics////////////////////////
        dg::direction diff_dir = m_p.diff_dir;
        dg::blas1::subroutine( [mu, tau, beta, diff_dir] DG_DEVICE (
                double N, double d0FN, double d1FN,
                          double d0BN, double d1BN,
                double U, double d0FU, double d1FU,
                          double d0BU, double d1BU,
                          double d0P, double d1P,
                double A, double d0A, double d1A,
                double b_2,
                double curvNabla0,  double curvNabla1,
                double curvKappa0,  double curvKappa1,
                double divCurvKappa,
                double& dtU
            )
            {
                dtU = 0;
                // velocity - upwind scheme
                double v0 = (b_1*d2P - b_2*d1P) + tau*curvNabla0 + mu*U*U*curvKappa0;
                double v1 = (b_2*d0P - b_0*d2P) + tau*curvNabla1 + mu*U*U*curvKappa1;
                double v2 = (b_0*d1P - b_1*d0P) + tau*curvNabla2 + mu*U*U*curvKappa2;
                double bp0 = 0., bp1 = 0., bp2 = 0.;
                if( beta != 0)
                {
                    bp0 = A * curvKappa0 + ( d1A*b_2 - d2A*b_1);
                    bp1 = A * curvKappa1 + ( d2A*b_0 - d0A*b_2);
                    bp2 = A * curvKappa2 + ( d0A*b_1 - d1A*b_0);

                    v0 += U * bp0;
                    v1 += U * bp1;
                    v2 += U * bp2;
                    //Q: doesn't U in U^2K_kappa and U b_perp create nonlinearity
                    //in velocity equation that may create shocks?
                    //A: we did some studies in the reconnection2d program and
                    //did not find shocks. LeVeque argues that for smooth
                    //solutions the upwind discretization should be fine but is
                    //wrong for shocks
                }
                // velocity - upwind scheme
                v0 += 2.*tau*curvKappa0;
                v1 += 2.*tau*curvKappa1;
                v2 += 2.*tau*curvKappa2;
                dtU += ( v0 > 0 ) ? -v0*d0BU : -v0*d0FU;
                dtU += ( v1 > 0 ) ? -v1*d1BU : -v1*d1FU;
                dtU += ( v2 > 0 ) ? -v2*d2BU : -v2*d2FU;

                double d0N = 0, d1N = 0, d2N = 0;
                if( diff_dir == dg::forward)
                    d0N = d0FN, d1N = d1FN, d2N = d2FN;
                else if( diff_dir == dg::backward)
                    d0N = d0BN, d1N = d1BN, d2N = d2BN;
                else
                    d0N = (d0FN+d0BN)/2., d1N = (d1FN+d1BN)/2.,
                        d2N = (d2FN+d2BN)/2.;
                // use centered derivatives
                double KappaN = curvKappa0*d0N+curvKappa1*d1N+curvKappa2*d2N;
                double KappaP = curvKappa0*d0P+curvKappa1*d1P+curvKappa2*d2P;

                dtU +=  - U * ( 2. * tau * KappaN / N + tau * divCurvKappa
                                + KappaP);
                if( beta != 0)
                {
                    double bpN = bp0 * d0N + bp1 * d1N + bp2 * d2N;
                    double bpP = bp0 * d0P + bp1 * d1P + bp2 * d2P;
                    dtU +=  - bpP/mu - tau/mu * bpN/N;
                }
                return;
            },
            //species depdendent
            density[i],   m_dFN[i][0], m_dFN[i][1],
                          m_dBN[i][0], m_dBN[i][1],
            velocity[i],  m_dFU[i][0], m_dFU[i][1],
                          m_dBU[i][0], m_dBU[i][1],
                          m_dP[i][0], m_dP[i][1],
            //aparallel
            apar, m_dA[0], m_dA[1],
            //magnetic parameters
            m_b[2],
            m_curvNabla[0], m_curvNabla[1],
            m_curvKappa[0], m_curvKappa[1],
            m_divCurvKappa, velocityDOT[i]
        );
    }
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix,
     Container>::compute_parallel_advection(
             const Container& velocityKM,
             const Container& velocityKP,
             const Container& densityM,
             const Container& density,
             const Container& densityP,
             Container& fluxM,
             Container& fluxP,
             std::string slope_limiter
             )
{
    // positive and negative velocities are defined wrt to the coordinate system
    // but we need it wrt to the b-field
    if( m_reversed_field)
    {
        dg::blas1::pointwiseDot( -1., velocityKM, m_vbm);
        dg::blas1::pointwiseDot( -1., velocityKP, m_vbp);
    }
    else
    {
        dg::blas1::copy( velocityKM, m_vbm);
        dg::blas1::copy( velocityKP, m_vbp);
    }
    dg::blas1::evaluate( fluxM, dg::equals(), dg::Upwind(),
            m_vbm, densityM, density);
    dg::blas1::evaluate( fluxP, dg::equals(), dg::Upwind(),
            m_vbp, density, densityP);
    if(slope_limiter != "none" )
    {
        // compute dn_k-1, dn_k, dn_k+1
        // By transforming dn to plus and minus planes
        dg::blas1::axpby( 1., densityP, -1., density, m_dNZ);

        m_fa( dg::geo::zeroForw, m_dNZ, m_dNP);
        m_fa( dg::geo::einsPlus, m_dNZ, m_dNPP);

        dg::blas1::axpby( 1., density, -1., densityM, m_dNZ);

        m_fa( dg::geo::zeroForw, m_dNZ, m_dNM);
        m_fa( dg::geo::einsMinus, m_dNZ, m_dNMM);

        // Let's keep the default boundaries of NEU
        // boundary values are (probably?) never used in the slope limiter branches
        //dg::blas1::copy(density, m_temp0); // save density
        //update_parallel_bc_2nd( m_fa, m_temp0, densityP, m_plus, dg::NEU, 0.);
        //dg::blas1::copy(density, m_temp0);
        //update_parallel_bc_2nd( m_fa, m_minus, densityM, m_temp0, dg::NEU, 0.);
        // dn is computed inside the limiter
        if( slope_limiter == "minmod")
        {
            dg::blas1::evaluate( fluxM, dg::plus_equals(),
                dg::SlopeLimiter<dg::MinMod>(), m_vbm,
                m_dNMM, m_dNM, m_dNP, 0.5, 0.5);
            dg::blas1::evaluate( fluxP, dg::plus_equals(),
                dg::SlopeLimiter<dg::MinMod>(), m_vbp,
                m_dNM, m_dNP, m_dNPP, 0.5, 0.5);
        }
        else if( slope_limiter == "vanLeer")
        {
            dg::blas1::evaluate( fluxM, dg::plus_equals(),
                dg::SlopeLimiter<dg::VanLeer>(), m_vbm,
                m_dNMM, m_dNM, m_dNP, 0.5, 0.5);
            dg::blas1::evaluate( fluxP, dg::plus_equals(),
                dg::SlopeLimiter<dg::VanLeer>(), m_vbp,
                m_dNM, m_dNP, m_dNPP, 0.5, 0.5);
        }
    }
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix,
     Container>::compute_parallel_flux(
             const Container& velocityKM,
             const Container& velocityKP,
             const Container& densityM,
             const Container& density,
             const Container& densityP,
             Container& fluxM,
             Container& fluxP,
             std::string slope_limiter
             )
{
    compute_parallel_advection( velocityKM, velocityKP, densityM, density, densityP,
            fluxM, fluxP, slope_limiter);
    dg::blas1::pointwiseDot( fluxM, velocityKM, fluxM);
    dg::blas1::pointwiseDot( fluxP, velocityKP, fluxP);
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix,
     Container>::compute_parallel_advection( const Container& velocity,
             const Container& minusST, const Container& plusST,
             Container& flux,
             std::string slope_limiter
             )
{
    if( m_reversed_field)
    {
        dg::blas1::pointwiseDot( -1., velocity, m_vbp);
    }
    else
    {
        dg::blas1::copy( velocity, m_vbp);
    }
    dg::blas1::evaluate( flux, dg::equals(), dg::Upwind(),
            m_vbp, minusST, plusST);
    if(slope_limiter != "none" )
    {
        // compute dn_k-1, dn_k, dn_k+1
        // By transforming dn to plus and minus planes
        dg::blas1::axpby( 1., plusST, -1., minusST, m_dN);
        m_fa( dg::geo::einsMinus, m_dN, m_dNM);
        m_fa( dg::geo::zeroForw,  m_dN, m_dNZ);
        m_fa( dg::geo::einsPlus,  m_dN, m_dNP);
        // Let's keep the default boundaries of NEU
        // boundary values are (probably?) never used in the slope limiter branches
        //update_parallel_bc_2nd( m_fa, m_minus, m_temp0, m_plus, dg::NEU, 0.);
        if( slope_limiter == "minmod")
        {
            dg::blas1::evaluate( flux, dg::plus_equals(),
                dg::SlopeLimiter<dg::MinMod>(), m_vbp,
                m_dNM, m_dNZ, m_dNP, 0.5, 0.5);
        }
        else if( slope_limiter == "vanLeer")
        {
            dg::blas1::evaluate( flux, dg::plus_equals(),
                dg::SlopeLimiter<dg::VanLeer>(), m_vbp,
                m_dNM, m_dNZ, m_dNP, 0.5, 0.5);
        }
    }
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix,
     Container>::compute_parallel_flux( const Container& velocity,
             const Container& minusST, const Container& plusST,
             Container& flux,
             std::string slope_limiter
             )
{
    compute_parallel_advection( velocity, minusST, plusST, flux, slope_limiter);
    dg::blas1::pointwiseDot( velocity, flux, flux);
}


template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::compute_parallel(
    std::array<std::array<Container,2>,2>& yp)
{
    for( unsigned i=0; i<2; i++)
    {
        // "velocity-staggered-fieldaligned"
        //// compute qhat
        //compute_parallel_flux( m_minusSTU[i], m_plusSTU[i],
        //        m_minusN[i], m_zeroN[i], m_plusN[i],
        //        m_minus, m_plus, m_p.slope_limiter);
        //// Now compute divNUb
        //dg::geo::ds_divCentered( m_faST, 1., m_minus, m_plus, 0.,
        //        m_divNUb[i]);
        //dg::blas1::axpby( -1., m_divNUb[i], 1., yp[0][i]);

        //// compute grad U2/2
        //dg::blas1::axpby( 0.25, m_minusU[i], 0.25, m_zeroU[i], m_minusSTU[i]);
        //dg::blas1::axpby( 0.25, m_zeroU[i],  0.25, m_plusU[i], m_plusSTU[i]);
        //compute_parallel_flux( m_minusSTU[i], m_plusSTU[i],
        //        m_minusU[i], m_zeroU[i], m_plusU[i],
        //        m_minus, m_plus,
        //        m_p.slope_limiter);
        //dg::geo::ds_centered( m_faST, -1., m_minus, m_plus, 1., yp[1][i]);
        //
        // "velocity-staggered"
        compute_parallel_flux( m_zeroU[i], m_minusSTN[i], m_plusSTN[i],
                m_temp0, m_p.slope_limiter);
        m_faST( dg::geo::zeroPlus,  m_temp0, m_plus);
        m_faST( dg::geo::einsMinus, m_temp0, m_minus);
        update_parallel_bc_1st( m_minus, m_plus, dg::NEU, 0.);
        dg::geo::ds_divCentered( m_faST, 1., m_minus, m_plus, 0., m_divNUb[i]);
        dg::blas1::axpby( -1., m_divNUb[i], 1., yp[0][i]);

        // compute fhat
        compute_parallel_flux( m_velocity[i], m_minusSTU[i], m_plusSTU[i],
                m_temp0, m_p.slope_limiter);
        m_faST( dg::geo::einsPlus, m_temp0, m_plus);
        m_faST( dg::geo::zeroMinus, m_temp0, m_minus);
        update_parallel_bc_1st( m_minus, m_plus, dg::NEU, 0.);
        dg::geo::ds_centered( m_faST, -0.5, m_minus, m_plus, 1, yp[1][i]);

        // Add density gradient and electric field
        double tau = m_p.tau[i], mu = m_p.mu[i], delta = m_fa.deltaPhi();
        dg::blas1::subroutine( [tau, mu, delta ]DG_DEVICE ( double& WDot,
                    double dsP, double QN, double PN, double bphi)
                {
                    WDot -= 1./mu*dsP;
                    WDot -= tau/mu*bphi*(PN-QN)/delta/2.*(1/PN + 1/QN);
                },
                yp[1][i], m_dsP[i], m_minusSTN[i], m_plusSTN[i], m_fa.bphi()
        );
    }
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::add_source_terms(
    std::array<std::array<Container,2>,2>& yp)
{
    if( m_source_rate != 0.0)
    {
        if( m_fixed_profile )
            dg::blas1::subroutine(
                [] DG_DEVICE ( double& result, double ne, double profne,
                    double source, double source_rate){
                    result = source_rate*source*(profne - ne);
                    },
                m_s[0][0], m_density[0], m_profne, m_source, m_source_rate);
        else
            dg::blas1::axpby( m_source_rate, m_source, 0., m_s[0][0]);
    }
    else
        dg::blas1::copy( 0., m_s[0][0]);
    // add prevention to get below lower limit
    if( m_minrate != 0.0)
    {
        // do not make lower forcing a velocity source
        // MW it may be that this form does not go well with the potential
        dg::blas1::transform( m_density[0], m_temp0, dg::PolynomialHeaviside(
                    m_minne-m_minalpha/2., m_minalpha/2., -1) );
        dg::blas1::transform( m_density[0], m_temp1, dg::PLUS<double>( -m_minne));
        dg::blas1::pointwiseDot( -m_minrate, m_temp1, m_temp0, 1., yp[0][0]);
        dg::blas1::transform( m_density[1], m_temp0, dg::PolynomialHeaviside(
                    m_minne-m_minalpha/2., m_minalpha/2., -1) );
        dg::blas1::transform( m_density[1], m_temp1, dg::PLUS<double>( -m_minne));
        dg::blas1::pointwiseDot( -m_minrate, m_temp1, m_temp0, 1., yp[0][1]);
    }

    //compute FLR corrections S_N = (1-0.5*mu*tau*Lap)*S_n
    dg::blas2::gemv( m_lapperpN, m_s[0][0], m_temp0);
    dg::blas1::axpby( 1., m_s[0][0], 0.5*m_p.tau[1]*m_p.mu[1], m_temp0, m_s[0][1]);
    // potential part of FLR correction S_N += -div*(mu S_n grad*Phi/B^2)
    dg::blas1::pointwiseDot( m_p.mu[1], m_s[0][0], m_binv, m_binv, 0., m_temp0);
    m_lapperpP.set_chi( m_temp0);
    m_lapperpP.symv( 1., m_potential[0], 1., m_s[0][1]);

    // S_U = - U S_N/N
    for(int i=0; i<2; i++)
    {
        // transform to adjoint plane and add to velocity source
        m_faST( dg::geo::zeroMinus, m_s[0][i], m_minus);
        m_faST( dg::geo::einsPlus,  m_s[0][i], m_plus);
        update_parallel_bc_1st( m_minus, m_plus, m_p.bcxN, 0.);
        dg::geo::ds_average( m_faST, 1., m_minus, m_plus, 0., m_temp0);
        dg::blas1::evaluate( m_s[1][i], dg::equals(), []DG_DEVICE(
                    double sn, double u, double n){ return -u*sn/n;},
                m_temp0, m_velocityST[i], m_densityST[i]);
    }
    //Add all to the right hand side
    dg::blas1::axpby( 1., m_s, 1.0, yp);
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::multiply_rhs_penalization(
        Container& yp)
{
    //mask right hand side in penalization region
    if( m_p.penalize_wall && m_p.penalize_sheath)
    {
        dg::blas1::subroutine( []DG_DEVICE(
            double& rhs, double wall, double sheath){
                rhs *= (1.0-wall-sheath);
            }, yp, m_wall, m_sheath);
    }
    else if( m_p.penalize_wall)
    {
        dg::blas1::subroutine( []DG_DEVICE( double& rhs, double wall){
                rhs *= (1.0-wall); }, yp, m_wall);
    }
    else if( m_p.penalize_sheath)
    {
        dg::blas1::subroutine( []DG_DEVICE( double& rhs, double sheath){
                rhs *= (1.0-sheath); }, yp, m_sheath);
    }
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::add_wall_and_sheath_terms(
        std::array<std::vector<Container>,4>& yp)
{
    // add sheath boundary conditions
    if( m_sheath_rate != 0)
    {
        ////density
        ////Here, we need to find out where "downstream" is
        //!! Simulations does not really work without
        for( unsigned i=0; i<2; i++)
        {
            //The coordinate automatically sees the reversed field
            //but m_plus and m_minus are defined wrt the angle coordinate
            if( m_reversed_field) //bphi negative (exchange + and -)
                dg::blas1::evaluate( m_temp0, dg::equals(), dg::Upwind(),
                     m_sheath_coordinate, m_plusN[i], m_minusN[i]);
            else
                dg::blas1::evaluate( m_temp0, dg::equals(), dg::Upwind(),
                     m_sheath_coordinate, m_minusN[i], m_plusN[i]);
            dg::blas1::pointwiseDot( m_sheath_rate, m_temp0, m_sheath, 1.,
                    yp[0][i]);
        }
        //compute sheath velocity
        if( "wall" == m_p.sheath_bc)
        {
            for( unsigned i=0; i<2; i++)
            {
                //dg::blas1::axpby( +m_sheath_rate*m_nwall, m_sheath, 1., yp[0][i] );
                dg::blas1::axpby( +m_sheath_rate*m_uwall, m_sheath, 1., yp[1][i] );
            }
        }
        else
        {
            //velocity c_s
            double cs = sqrt(1.+m_p.tau[1]), sheath_rate = m_sheath_rate;
            if( "insulating" == m_p.sheath_bc)
            {
                // u_e,sh = s*sqrt(1+tau) Ni/ne
                dg::blas1::evaluate( yp[1][0], dg::plus_equals(),
                        [cs, sheath_rate]DG_DEVICE( double sheath_coord, double
                            sheath, double ne, double ni) {
                            return cs*sheath_rate*sheath_coord*ni/ne*sheath;
                        },
                        m_sheath_coordinate, m_sheath, m_densityST[0],
                        m_densityST[1]);
            }
            else // "bohm" == m_p.sheath_bc
            {
                //u_e,sh = s*1/sqrt(|mu_e|2pi) exp(-phi)
                double mue = fabs(m_p.mu[0]), tau = m_p.tau[1];
                dg::blas1::evaluate( yp[1][0], dg::plus_equals(),
                    [mue, sheath_rate, tau]DG_DEVICE(
                        double sheath_coord, double sheath, double phi) {
                        return sheath_rate * sheath_coord * sheath *
                            sqrt(1.+tau) * exp(-phi) / sqrt( mue*2.*M_PI);
                    },
                    m_sheath_coordinate, m_sheath, m_potentialST[0]);
            }
            // u_i,sh = s*sqrt(1+tau)
            dg::blas1::pointwiseDot( sheath_rate*cs,
                    m_sheath, m_sheath_coordinate, 1.,  yp[1][1]);
        }
    }
    // add wall boundary conditions
    if( m_wall_rate != 0)
    {
        for( unsigned i=0; i<2; i++)
        {
            dg::blas1::axpby( +m_wall_rate*m_nwall, m_wall, 1., yp[0][i] );
            dg::blas1::axpby( +m_wall_rate*m_uwall, m_wall, 1., yp[1][i] );
        }
    }
}

#ifndef WITH_NAVIER_STOKES
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::operator()(
    double t,
    const std::array<std::vector<Container>,4>& y,
    std::array<std::vector<Container>,4>& yp)
{
    m_called++;
#ifdef MPI_VERSION
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
    //DG_RANK0 std::cout << "## time "<<time<<" dt "<<dt<<" t_out "<<t_output<<" step "<<step<<" failed "<<var.nfailed<<"\n";
    DG_RANK0 std::cout << "## time "<<t<<"\n";

    dg::Timer timer;
    double accu = 0.;//accumulated time
    timer.tic();


    std::vector<Container>&  density    = y[0];
    std::vector<Container>&  wST        = y[1];
    std::vector<Container>&  pperp      = y[2];
    std::vector<Container>&  pparaST    = y[3];

    // First compute potentials phi and apar
#if FELTORPERP == 1

    // set m_potential[0]
    m_solvers.compute_phi( t, density, pperp, multiply_rhs_penalization, m_phi);

#else


#endif
    timer.toc();
    accu += timer.diff();
    DG_RANK0 std::cout << "## Compute phi                       took "
                       << timer.diff()<<"s\t A: "<<accu<<"s\n";
    timer.tic( );
    //Compute m_densityST, m_minusSTN, m_plusSTN
    update_staggered_density( t, m_density);

    // Compute m_aparST and m_velocityST if necessary
    if( m_p.beta != 0)
    {
        m_solvers.compute_aparST( t, m_densityST, wST, m_aparST, true);
    }
    // Compute apar
    m_faST( dg::geo::einsMinus, m_aparST, m_minus);
    m_faST( dg::geo::zeroPlus,  m_aparST, m_plus);
    update_parallel_bc_1st( m_minus, m_plus, m_p.bcxA, 0.);
    dg::blas1::axpby( 0.5, m_minus, 0.5, m_plus, m_apar);
    m_old_apar.update( t, m_apar);

    // main species loop
    for( unsigned s=0; s<m_p.num_species; s++)
    {

        m_upToDate[s] = false;
        //Compute m_velocityST, m_minusSTU, m_plusST, m_minusU, m_zeroU, m_plusU
        update_velocity( t, wST, m_aparST, m_velocityST, s);

        // Compute the three potentials
        m_solvers.compute_psi( t, 

        //Compute m_minusN, m_zeroN, m_plusN
        update_staggered_density_and_phi( t, m_density, m_potential);

        timer.toc();
        accu += timer.diff();
        DG_RANK0 std::cout << "## Compute phi and psi ST            took "
                           << timer.diff()<<"s\t A: "<<accu<<"s\n";
        timer.tic( );


#if FELTORPERP == 1

        // Set perpendicular dynamics in yp
        compute_perp_density(  t, m_density, m_velocity, m_potential, m_apar,
                yp[0]);
        compute_perp_velocity( t, m_densityST, m_velocityST, m_potentialST,
                m_aparST, yp[1]);
        compute_perp_pperp( t, m_densityST, m_velocityST, m_potentialST,
                m_aparST, yp[1]);
        compute_perp_para( t, m_densityST, m_velocityST, m_potentialST,
                m_aparST, yp[1]);

#else

        dg::blas1::copy( 0., yp);

#endif

        timer.toc();
        accu += timer.diff();
        DG_RANK0 std::cout << "## Compute perp dynamics             took "
                           << timer.diff() << "s\t A: "<<accu<<"s\n";
        timer.tic();

        // Add parallel dynamics
#if FELTORPARALLEL == 1

        compute_parallel( yp);

#endif
#if FELTORPERP == 1
        //------------------Add Resistivity--------------------------//
        double eta = m_p.eta, mu0 = m_p.mu[0], mu1 = m_p.mu[1];
        dg::blas1::subroutine( [eta,mu0,mu1] DG_DEVICE (
                double ne, double ni,
                double ue, double ui, double& dtUe, double& dtUi){
                    double current = ni*ui-ne*ue;
                    dtUe += -eta/mu0 * current;
                    dtUi += -eta/mu1 * ne/ni * current;
                },
            m_densityST[0], m_densityST[1],
            m_velocityST[0], m_velocityST[1], yp[1][0], yp[1][1]);
#endif
    } // species loops

    if( !m_p.partitioned)
    {
        // explicit and implicit timestepper
        add_implicit_density( t, m_density, 1., yp[0]);
        add_implicit_velocityST( t, m_densityST, m_velocityST, 1., yp[1]);
        add_implicit_perp( t, m_perp, 1., yp[2]);
        add_implicit_pparaST( t, m_pparaST, m_velocityST, 1., yp[3]);
    }
    else
    {
        // partitioned means imex timestepper
        for( unsigned i=0; i<4; i++)
        {
            for( unsigned s=0; s<m_p.num_species; s++)
                multiply_rhs_penalization( yp[i][s]); // F*(1-chi_w-chi_s)
        }
    }

    add_wall_and_sheath_terms( yp);
    //Add source terms
    // set m_s
    add_source_terms( yp );

    timer.toc();
    accu += timer.diff();
    #ifdef MPI_VERSION
        if(rank==0)
    #endif
    std::cout << "## Add parallel dynamics and sources took "<<timer.diff()
              << "s\t A: "<<accu<<"\n";
}
template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::add_implicit_density(
    double t,
    const std::vector<Container>& density,
    double beta,
    std::vector<Container>& yp)
{
#if FELTORPERP == 1
    for( unsigned i=0; i<m_p.num_species; i++)
        compute_perp_diffusiveN( 1., density[i], m_temp0,
                m_temp1, beta, yp[i]);
#else
    dg::blas1::scal( yp, beta);
#endif
#if FELTORPARALLEL == 1
    for( unsigned i=0; i<m_p.num_species; i++)
    {
        if( m_p.nu_parallel_n > 0)
        {
            dg::geo::dssd_centered( m_fa, m_p.nu_parallel_n,
                    m_minusN[i], m_zeroN[i], m_plusN[i], 1., yp[i]);
        }
    }
#endif
    for( unsigned i=0; i<m_p.num_species; i++)
    {
        multiply_rhs_penalization( yp[i]); // F*(1-chi_w-chi_s)
        dg::blas1::pointwiseDot( -m_wall_rate, m_wall, density[i],
            -m_sheath_rate, m_sheath, density[i], 1., yp[i]); // -r N
    }
}

template<class Geometry, class IMatrix, class Matrix, class Container>
template<size_t N>
void Explicit<Geometry, IMatrix, Matrix, Container>::add_implicit_velocityST(
    double t,
    const std::vector<Container>& densityST,
    const std::vector<Container>& velocityST,
    double beta,
    std::vector<Container>& yp)
{
    // velocityST[s] := u_s^dagger
#if FELTORPERP == 1
    for( unsigned i=0; i<m_p.num_species; i++)
    {
        compute_perp_diffusiveU( 1., velocityST[i], densityST[i], m_temp0,
                m_temp1, m_dFU[i][0], m_dFU[i][1], beta, yp[i]);
    }
#else
    dg::blas1::scal( yp, beta);
#endif
#if FELTORPARALLEL == 1
    for( unsigned i=0; i<m_p.num_species; i++)
    {
        // Add parallel viscosity
        if( m_p.nu_parallel_u[i] > 0)
        {
            dg::geo::dssd_centered( m_fa, m_p.nu_parallel_u[i],
                    m_minusU[i], m_zeroU[i], m_plusU[i], 0., m_temp0);
            if( !m_modify_diff)
                dg::blas1::pointwiseDivide( 1., m_temp0, densityST[i], 1., yp[i]);
            else
                dg::blas1::axpby( 1., m_temp0, 1., yp[i]);
        }
        double nu = m_p.nu_parallel_n;
        if( m_modify_diff)
            nu += m_p.nu_parallel_u[i];
        if( nu > 0)
        {
            // Add density gradient correction
            double delta = m_fa.deltaPhi();
            dg::blas1::subroutine( [delta, nu]DG_DEVICE ( double& WDot,
                        double QN, double PN, double UM, double U0, double UP,
                        double bphi)
                    {
                        //upwind scheme
                        double nST = (PN+QN)/2.;
                        double current = -nu*bphi*(PN-QN)/delta/nST;
                        if( current > 0)
                            WDot += - current*bphi*(U0-UM)/delta;
                        else
                            WDot += - current*bphi*(UP-U0)/delta;

                    },
                    yp[i], m_minusSTN[i], m_plusSTN[i], m_minusU[i], m_zeroU[i],
                    m_plusU[i], m_fa.bphi()
            );
        }
    }
#endif
    for( unsigned i=0; i<m_p.num_species; i++)
    {
        multiply_rhs_penalization( yp[i]); // F*(1-chi_w-chi_s)
        dg::blas1::pointwiseDot( -m_wall_rate, m_wall, velocityST[i],
            -m_sheath_rate, m_sheath, velocityST[i], 1., yp[i]); // -r U
    }
}

//#else // WITH_NAVIER_STOKES
//#include "../navier_stokes/navier_stokes.h"
//#endif // WITH_NAVIER_STOKES

} //namespace thermal
