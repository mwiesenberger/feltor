#pragma once

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "parameters.h"
#include "solvers.h"
#include "../feltor/common.h"
#include "perpendicular.h"
#include "parallel.h"

namespace thermal
{

template< class Geometry, class IMatrix, class Matrix, class Container >
struct Explicit
{
    // y[0] == density
    // y[1] == pperp
    // y[2] == ppara
    // y[3] == wST
    // y[4] == qperpST
    // y[5] == qparaST
    using vector = std::array<std::vector<Container>,4>;
    using container = Container;
    Explicit( const Geometry& g, thermal::Parameters p,
        dg::geo::TokamakMagneticField mag, dg::file::WrappedJsonValue js );

    void operator()( double t,
        const std::array<std::vector<Container>,6>& y,
        std::array<std::vector<Container>,6>& yp);
    /// ///////////////////DIAGNOSTIC MEMBERS //////////////////////
    /// 3d output
    const Geometry& grid() const {
        return m_solvers.grid();
    }
    const Container& potential() const {
        return m_phi;
    }
    const Container& aparallel() const {
        return m_apar;
    }
    /// 2d static output
    const std::array<Container, 2> & curvNabla () const {
        return m_perp.curvNabla();
    }
    const std::array<Container, 2> & curvKappa () const {
        return m_perp.curvKappa();
    }
    const Container& divCurvKappa() const {
        return m_perp.divCurvKappa();
    }
    const Container& bphi( ) const { return m_perp.bphi(); }
    const Container& divb( ) const { return m_perp.divb(); }
    const Container& binv( ) const { return m_perp.binv(); }
    const Container& bhatgB( ) const { return m_perp.bhatgB(); } // \pm 1/B
    const Container& get_wall() const{
        return m_wall;
    }
    const Container& get_sheath() const{
        return m_parallel.get_sheath();
    }
    const Container& get_sheath_coordinate() const{
        return m_parallel.get_sheath_coordinate();
    }

    // Chanes after each species
    const std::array<Container, 4>& psi() const {
        return m_psi;
    }



    const Container& uE2() const {
        return m_UE2;
    }
    const Container&  gammaNi() const{
        return m_old_gammaN.head();
    }
    const Container&  gammaPhi() const{
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
        dg::blas1::transform( density, temp0, dg::PLUS<double>(-m_p.nbc[s]));
        dg::blas2::symv( alpha, m_lapperpN, temp0, beta, result);
    }
    void compute_lapMperpU (int i, Container& result)
    {
        dg::blas2::symv( m_lapperpU, m_velocity[i], result);
    }
    void compute_lapMperpP (int i, Container& result)
    {
        m_lapperpP.set_chi( 1.);
        dg::blas2::gemv( m_lapperpP, m_psi[i], result);
    }
    void compute_lapMperpA ( Container& result)
    {
        // only if lapperpU has same direction as lapperpP
        dg::blas2::gemv( m_lapperpU, m_apar, result);
    }
    const std::array<Container,2> get_bperp( )
    {
        update_diag();
        return m_bperp;
    }
    const Container& get_source() const{
        return m_source;
    }
    const Container& get_source_prof() const{
        return m_profne;
    }
    void compute_perp_diffusiveN( double alpha, const Container& density,
            Container& temp0, Container& temp1, double beta, Container& result )
    {
        // density = full N
        // result = alpha Lambda_N + beta result
        if( m_p.nu_perp_n > 0)
        {
            dg::blas1::transform( density, temp0, dg::PLUS<double>(-m_p.nbc[s]));
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
            dg::blas1::pointwiseDivide( -alpha*m_p.nu_perp_u, temp0, density, beta, result);
        }
        else
            dg::blas1::scal( result, beta);
        double nu = m_p.nu_perp_n;
        if( nu > 0 )
        {

            dg::blas1::transform( density, temp0, dg::PLUS<double>(-m_p.nbc[s]));
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

    void compute_parallel_diffusiveN( int i, Container& result)
    {
        dg::blas1::axpby( m_p.nu_parallel_n, lapParN(i), 0., result);
    }
    void compute_parallel_diffusiveU( int i, Container& result)
    {
        double nu = m_p.nu_parallel_n;
        if( nu > 0)
        {
            dg::blas1::pointwiseDot( dsN(i), dsU(i), result);
            dg::blas1::pointwiseDivide( nu, result, density(1), 0., result);
        }
        else
            dg::blas1::copy( 0, result);
        if( m_p.nu_parallel_u[i] > 0)
        {
            dg::blas1::pointwiseDivide( m_p.nu_parallel_u[i], lapParU(i), density(i), 1., result);
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
        dg::blas2::symv( -alpha, m_multi_pol[0], m_psi[0], beta, result);
    }
    void compute_source_pol( double alpha, const Container& density, Container& temp, double beta, Container& result)
    {
        // we don't want jumps in phi in here so we use lapperpP
        dg::blas1::pointwiseDot( m_p.mu[1], density, m_binv, m_binv, 0., temp);
        m_lapperpP.set_chi( temp);
        dg::blas2::symv( -alpha, m_lapperpP, m_psi[0], beta, result);
    }
    unsigned called() const { return m_called;}

    /// //////////////////////DIAGNOSTICS END////////////////////////////////
    void update_diag(){
        // assume m_density, m_potential, m_velocity, m_velocityST, m_apar
        // compute dsN, dsU, dsP, lapParU, dssU and perp derivatives
        if( !m_upToDate)
        {
            // update m_dN, m_dU, m_dP, m_dA
            update_perp_derivatives( m_density, m_velocity, m_psi, m_apar);
            for( unsigned i=0; i<2; i++)
            {
                // density m_dsN, m_lapParN
                m_fa( dg::geo::einsMinus, m_density[i], m_minus);
                m_fa( dg::geo::zeroForw,  m_density[i], m_zero);
                m_fa( dg::geo::einsPlus,  m_density[i], m_plus);
                update_parallel_bc_2nd( m_fa, m_minus, m_zero, m_plus,
                        m_p.bcxN, m_p.bcxN == dg::DIR ? m_p.nbc[s] : 0.);
                dg::geo::ds_centered( m_fa, 1., m_minus, m_plus, 0., m_dsN[i]);
                dg::geo::dssd_centered( m_fa, 1.,
                        m_minus, m_zero, m_plus, 0., m_lapParN[i]);
                // potential m_dsP
                m_fa( dg::geo::einsMinus, m_psi[i], m_minus);
                m_fa( dg::geo::einsPlus,  m_psi[i], m_plus);
                update_parallel_bc_2nd( m_fa, m_minus, m_psi[i], m_plus,
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
    void set_source(
        bool fixed_profile, // cannot be mixed among species or equations due to quasineutrality and transformation
        const std::array<std::vector<double>,3>& source_rate,
        const std::array<std::vector<dg::x::HVec>,3>& profile, // for influx this can be ignored
        const std::array<std::vector<dg::x::HVec>,3>& source,  // for fixed profile this contains damping
        const std::vector<double>& minne,
        double minrate,
        const std::vector<double>& minalpha)
    {
        m_fixed_profile = fixed_profile;
        m_source_rate = source_rate;
        for( unsigned u=0; u<3; u++)
        {
            m_profile.resize( m_p.num_species);
            m_source.resize( m_p.num_species);
            for( unsigned s=0; s<m_p.num_species; s++)
            {
                m_profile[u][s] = profile[u][s];
                m_source[u][s]  = source[u][s];
            }
        }
        m_minne = minne;
        m_minrate = minrate;
        m_minalpha = minalpha;
    }
    void set_wall(const Container& wall)
    {
        dg::blas1::copy( wall, m_wall);
    }
    void set_sheath(double sheath_rate, const Container& sheath,
            const Container& sheath_coordinate)
    {
        m_parallel.set_sheath( sheath_rate, sheath, sheath_coordinate);
    }
    void add_wall_terms(   std::array<std::vector<Container>,6>& yp);
    void add_source_terms( std::array<std::vector<Container>,6>& yp);
    const dg::geo::Fieldaligned<Geometry, IMatrix, Container>& fieldaligned() const
    {
        return m_parallel.fieldaligned();
    }
    // Called in init
    void transform_density_pperp( double mus, double zs, const Container& density, const Container& pperp, const Container& phi,
        Container& gydensity, Container& gypperp) // can be called inplace!
    {
        m_perp.transform_density_pperp( mus, zs, density, pperp, phi, gydensity, gypperp);
    }
  private:
    PerpDynamics<Geometry, IMatrix, Matrix, Container> m_perp;
    ParallelDynamics<Geometry, IMatrix, Matrix, Container> m_parallel;
    ThermalSolvers<Geometry, Matrix, Container> m_solvers;
    Collisions<Geometry, IMatrix, Matrix, Container> m_collisions;

    //
    Container m_phi, m_apar, m_aparST; // same for all species

    Container m_temp0, m_temp1;
    Container m_wall;
    dg::Extrapolation<Container> m_old_apar;
    std::array<Container,4> m_psi; // different for each species

    std::array<std::vector<Container>,3> m_s; // source terms for all species
    std::array<std::vector<Container>,3> m_source, m_profile; // (physical) source terms for all species

    const thermal::Parameters m_p;
    const dg::file::WrappedJsonValue m_js;
    std::vector<bool> m_upToDate;

    std::array<std::vector<double>,3> m_source_rate;
    bool m_nixed_profile = false;

    std::vector<double> m_minne, m_minalpha;
    double m_minrate  = 0.;
    unsigned m_called = 0;

};

template<class Grid, class IMatrix, class Matrix, class Container>
Explicit<Grid, IMatrix, Matrix, Container>::Explicit( const Grid& g,
    thermal::Parameters p, dg::geo::TokamakMagneticField mag,
    dg::file::WrappedJsonValue js
    ):
    m_perp( g, p, mag, js),
    m_parallel( g, p, mag, js),
    m_solvers( g, p, mag, js),
    m_collisions( g, p, mag, js),
    m_old_apar( 2, dg::evaluate( dg::zero, g)),
    m_p(p), m_js(js)
{
    //--------------------------init vectors to 0-----------------//
    dg::assign( dg::evaluate( dg::zero, g), m_phi );
    m_apar = m_aparST = m_phi;
    m_temp0 = m_temp1 = m_phi;
    for( int i=0; i<4; i++)
        m_psi[i] = m_phi;

    //--------------------------Construct-------------------------//
    m_upToDate.resize( m_p.num_species);
    for( unsigned s=0; s<m_p.num_species; s++)
        m_upToDate[s] = false;

}



template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::add_source_terms(
    const std::array<std::vector<Container>,6>& y,
    std::array<std::vector<Container>,6>& yp)
{
    for( unsigned s=0; s<m_p.num_species; s++)
    {
        // First, add minimum density source
        // add prevention to get below lower limit
        if( m_minrate != 0.0)
        {
            // do not make lower forcing a velocity source
            // MW it may be that this form does not go well with the potential
            dg::blas1::transform( y[0][s], m_temp0, dg::PolynomialHeaviside(
                        m_minne[s]-m_minalpha[s]/2., m_minalpha[s]/2., -1) );
            dg::blas1::transform( y[0][s], m_temp1, dg::PLUS<double>( -m_minne[s]));
            dg::blas1::pointwiseDot( -m_minrate, m_temp1, m_temp0, 1., yp[0][s]);
        }

        if( m_fixed_profile )
        {
            // If fixed profile transform profiles first
            transform_density_pperp( m_p.mu[s], m_p.z[s], m_profile[0][s],  m_profile[1][s], m_phi,
                m_s[0][s], m_s[1][s]);
            dg::blas1::copy( m_profile[2][s], m_s[2][s]);
            for( unsigned u=0; u<3; u++)
            {
                if( m_source_rate[u][s] != 0.0)
                {
                    dg::blas1::subroutine(
                        [] DG_DEVICE ( double& result, double ne, double profne,
                            double source, double source_rate){
                            result = source_rate*source*(profne - ne);
                            },
                        m_s[u][s], y[u][s], m_s[u][s], m_source[u][s], m_source_rate[u][s]);
                }
                else
                    dg::blas1::copy( 0., m_s[u][s]);
            }
        }
        else // influx
        {
            for( unsigned u=0; u<3; u++)
            {
                if( m_source_rate[u][s] != 0.0)
                {
                    dg::blas1::axpby( m_source_rate[u][s], m_source[u][s], 0., m_s[u][s]);
                }
                else
                    dg::blas1::copy( 0., m_s[u][s]);
            }
            // if influx transform sources last
            transform_density_pperp( m_p.mu[s], m_p.z[s], m_s[0][s],  m_s[1][s], m_phi,
                m_s[0][s], m_s[1][s]);
        }
    }
    //Add all to the right hand side
    for( unsigned u=0; u<3; u++)
        for( unsigned s=0; s<m_p.num_species; s++)
            dg::blas1::axpby( 1., m_s[u][s], 1.0, yp[u][s]);
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::add_wall_terms(
        std::array<std::vector<Container>,6>& yp)
{
    // add wall boundary conditions
    if( m_p.wall_rate != 0)
    {
        for( unsigned s=0; s<m_p.num_species; s++)
        {
            dg::blas1::axpby( +m_p.wall_rate*m_p.nwall[s], m_wall, 1., yp[0][s] );
            dg::blas1::axpby( +m_p.wall_rate*m_p.twall*m_p.nwall[s], m_wall, 1., yp[1][s] );
            dg::blas1::axpby( +m_p.wall_rate*m_p.twall*m_p.nwall[s], m_wall, 1., yp[2][s] );
            dg::blas1::axpby( +m_p.wall_rate*m_p.uwall, m_wall, 1., yp[3][s] );
            dg::blas1::axpby( +m_p.wall_rate*m_p.qwall, m_wall, 1., yp[4][s] );
            dg::blas1::axpby( +m_p.wall_rate*m_p.qwall, m_wall, 1., yp[5][s] );
        }
    }
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Explicit<Geometry, IMatrix, Matrix, Container>::operator()(
    double t,
    const std::array<std::vector<Container>,6>& y,
    std::array<std::vector<Container>,6>& yp)
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
    std::vector<Container>&  pperp      = y[1];
    std::vector<Container>&  ppara      = y[2];
    std::vector<Container>&  wST        = y[3];
    std::vector<Container>&  qperpST    = y[4];
    std::vector<Container>&  qparaST    = y[5];

    // First compute potentials phi and apar

    m_solvers.compute_phi( t, density, pperp, m_phi, m_p.penalize_wall,
        m_wall, m_p.penalize_sheath, m_sheath);

    timer.toc();
    accu += timer.diff();
    DG_RANK0 std::cout << "## Compute phi                       took "
                       << timer.diff()<<"s\t A: "<<accu<<"s\n";
    timer.tic( );
    //Compute m_densityST, m_minusSTN, m_plusSTN for all species
    m_parallel.update_staggered_density( t, density);

    // Compute m_aparST and m_velocityST if necessary
    if( m_p.beta != 0)
    {
        m_solvers.compute_aparST( t, m_parallel.get_staggered_density(), wST,
            m_aparST, true);
    }
    m_parallel.update_apar( m_aparST, m_apar);
    // Compute apar, bperp, and divbperp
    m_old_apar.update( t, m_apar);

    // main species loop
    for( unsigned s=0; s<m_p.num_species; s++)
    {
        m_upToDate[s] = false;

        // Compute the three potentials
        m_solvers.compute_psi( t, density, pperp, m_phi, m_psi, s);
        m_parallel.update_quantities( s, m_aparST, m_psi, y);

        // Set perpendicular dynamics in yp
        m_perp.add_perp_densities_advection(  s, m_apar, m_psi, y,
            m_parallel.get_q(), yp);
        m_perp.add_perp_velocities_advection( s, m_aparST,
            m_parallel.get_psiST(), y, m_parallel.get_qST(), yp);


        timer.toc();
        accu += timer.diff();
        DG_RANK0 std::cout << "## Compute perp dynamics             took "
                           << timer.diff() << "s\t A: "<<accu<<"s\n";
        timer.tic();

        // Add parallel dynamics
        m_parallel.add_para_densities_advection( s, m_psi, y, yp);
        m_parallel.add_para_velocities_advection( s, y, yp);

        m_parallel.add_perp_densities_diffusion( s, y, m_parallel.get_q(), yp);
        m_parallel.add_perp_velocities_diffusion( s, y, m_parallel.get_qST(), yp);
        m_parallel.add_para_densities_diffusion( s, y, yp);
        m_parallel.add_para_velocities_diffusion( s, y, yp);

        m_parallel.add_sheath_neumann_terms( s, yp);

        m_collisions.gather_quantities( s, m_parallel.get_q(), m_parallel.get_qST() )

        timer.toc();
        accu += timer.diff();
        DG_RANK0 std::cout << "## Compute para dynamics             took "
                           << timer.diff() << "s\t A: "<<accu<<"s\n";
        timer.tic();

    } // species loops

    m_collisions.add_collisions( m_parallel.get_qST()[0], m_aparST, y, yp);
    m_parallel.add_sheath_velocity_terms( m_collisions.get_tparaST(), yp);

    for( unsigned u=0; u<6; u++)
    {
        for( unsigned s=0; s<m_p.num_species; s++)
        {
            common::multiply_rhs_penalization( yp[u][s], m_p.penalize_wall,
                m_wall, m_p.penalize_sheath, m_sheath); // F*(1-chi_w-chi_s)
            if( u == 3)
                // Apply to U, not W
                dg::blas1::pointwiseDot( -m_wall_rate, m_wall, m_parallel.get_qST[3],
                    -m_sheath_rate, m_sheath, m_parallel.get_qST[3], 1., yp[u][s]);
            else
                dg::blas1::pointwiseDot( -m_wall_rate, m_wall, y[u][s],
                    -m_sheath_rate, m_sheath, y[u][s], 1., yp[u][s]);
        }
    }

    add_wall_terms( yp);
    //Add source terms
    // set m_s
    add_source_terms( yp );
    m_parallel.add_velocity_source_term( m_s[0], wST, m_aparST, yp);

    timer.toc();
    accu += timer.diff();
    DG_RANK0 std::cout << "## Add parallel dynamics and sources took "<<timer.diff()
              << "s\t A: "<<accu<<"\n";
}



} //namespace thermal
