#pragma once

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"

#include "../feltor/common.h"

#include "parameters.h"
#include "solvers.h"
#include "perpendicular.h"
#include "parallel.h"
#include "collisions.h"
#include "sources.h"


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
        return m_sources.get_wall();
    }
    const Container& get_sheath() const{
        return m_parallel.get_sheath();
    }
    const Container& get_sheath_coordinate() const{
        return m_parallel.get_sheath_coordinate();
    }

    unsigned called() const { return m_called;}
    const Container& get_source(unsigned u, unsigned s) const{
        return m_sources.get_source( u,s);
    }
    const Container& get_source_prof(unsigned u, unsigned s) const{
        return m_sources.get_source_prof( u,s);
    }

    /*
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
        return m_ss[0][i];
    }
    const Container& velocity(int i)const{
        return m_velocity[i];
    }
    const Container& velocity_source(int i){
        update_diag();
        return m_ss[1][i];
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
        dg::blas2::symv( m_dxF_N, m_ss[0][i], gradS[0]);
        dg::blas2::symv( m_dyF_N, m_ss[0][i], gradS[1]);
    }
    void compute_dot_aparallel( Container& tmp) const {
        m_old_apar.derive( tmp);
    }

    //volume with dG weights
    const Container& vol3d() const { return m_lapperpN.weights();}
    const Container& weights() const { return m_lapperpN.weights();}
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
                dg::blas1::evaluate( m_ss[1][i], dg::equals(), []DG_DEVICE(
                            double sn, double u, double n){ return -u*sn/n;},
                        m_ss[0][i], m_velocity[i], m_density[i]);
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
    */
    void set_source(
        bool fixed_profile, // cannot be mixed among species or equations due to quasineutrality and transformation
        const std::array<std::vector<double>,3>& source_rate,
        const std::array<std::vector<dg::x::HVec>,3>& profile, // for influx this can be ignored
        const std::array<std::vector<dg::x::HVec>,3>& source,  // for fixed profile this contains damping
        const std::vector<double>& minne,
        double minrate,
        const std::vector<double>& minalpha)
    {
        m_sources.set_source( fixed_profile, source_rate, profile, source,
            minne, minrate, minalpha);
    }
    void set_wall(const Container& wall)
    {
        m_sources.set_wall( wall);
    }

    void set_sheath(double sheath_rate, const Container& sheath,
            const Container& sheath_coordinate)
    {
        m_parallel.set_sheath( sheath_rate, sheath, sheath_coordinate);
    }
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
    Sources<Geometry, IMatrix, Matrix, Container> m_sources;

    dg::Extrapolation<Container> m_old_apar; // for diagnostics

    std::map<std::string, std::vector<Container>> m_q;

    const std::vector<std::string> q_names = {
        "N", "N 0", "N +1", "N -1",
        "Tperp", "Tperp 0", "Tperp +1", "Tperp -1",
        "Tpara", "Tpara 0", "Tpara +1", "Tpara -1",
        "Psi0",
        "Psi1", "ds Psi1",
        "Psi2",
        "Psi3",
        "U", "U +1/2", "U -1/2",
        "Uperp",
        "Upara",
        "Qperp +1/2", "Qperp -1/2",
        "Qpara +1/2", "Qpara -1/2",
        // Staggered variables
        "ST N", "ST N +1/2", "ST N -1/2",
        "ST Tperp", "ST ds Tperp",
        "ST Tpara", "ST ds Tpara",
        "ST Pperp +1/2", "ST Pperp -1/2",
        "ST Ppara +1/2", "ST Ppara -1/2",
        "ST Psi0", "ST ds Psi0",
        "ST Psi1", "ST ds Psi1",
        "ST Psi2",
        "ST Psi3",
        "ST U", "ST U 0", "ST U +1", "ST U -1",
        "ST Uperp", "ST Uperp 0", "ST Uperp +1", "ST Uperp -1",
        "ST Upara", "ST Upara 0", "ST Upara +1", "ST Upara -1",
        // perp derivatives
        "dxF N", "dxB N", "dyF N", "dyB N",
        "dxF Pperp", "dxB Pperp", "dyF Pperp", "dyB Pperp",
        "dxF Ppara", "dxB Ppara", "dyF Ppara", "dyB Ppara",
        "dx Tperp", "dy Tperp",
        "dx Tpara", "dy Tpara",
        "dx Psi0", "dy Psi0",
        "dx Psi1", "dy Psi1",
        "dx Psi2", "dy Psi2",
        "dx U", "dy U",
        "dx Uperp", "dy Uperp",
        "dx Upara", "dy Upara",
        // Staggered perp derivatives
        "ST dx N", "ST dy N",
        "ST dx Tperp", "ST dy Tperp",
        "ST dx Tpara", "ST dy Tpara",
        "ST dx Psi0", "ST dy Psi0",
        "ST dx Psi1", "ST dy Psi1",
        "ST dx Psi2", "ST dy Psi2",
        "ST dxF U", "ST dxB U", "ST dyF U", "ST dyB U",
        "ST dxF Qperp", "ST dxB Qperp", "ST dyF Qperp", "ST dyB Qperp",
        "ST dxF Qpara", "ST dxB Qpara", "ST dyF Qpara", "ST dyB Qpara"
    };

    Container m_phi, m_apar, m_dxapar, m_dyapar, m_aparST, m_dxaparST, m_dyaparST; // same for all species

    const thermal::Parameters m_p;
    const dg::file::WrappedJsonValue m_js;
    std::vector<bool> m_upToDate;

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
    m_sources( g, p, mag, js),
    m_old_apar( 2, dg::evaluate( dg::zero, g)),
    m_p(p), m_js(js)
{
    //--------------------------init vectors to 0-----------------//
    dg::assign( dg::evaluate( dg::zero, g), m_phi );
    m_apar = m_dxapar = m_dyapar = m_dxaparST = m_dyaparST = m_aparST = m_phi;
    m_q["N"] = std::vector<Container>( m_p.num_species, m_phi);
    for( auto name : q_names)
        m_q[name] = m_q["N"];
    //--------------------------Construct-------------------------//
    m_upToDate.resize( m_p.num_species);
    for( unsigned s=0; s<m_p.num_species; s++)
        m_upToDate[s] = false;

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


    const std::vector<Container>&  density    = y[0];
    const std::vector<Container>&  pperp      = y[1];
    const std::vector<Container>&  ppara      = y[2];
    const std::vector<Container>&  wST        = y[3];
    const std::vector<Container>&  qperpST    = y[4];
    const std::vector<Container>&  qparaST    = y[5];


    // 1. Transform Pperp, Ppara to Tperp, Tpara for all species
    dg::blas1::copy( y[0], m_q["N"]);
    dg::blas1::pointwiseDivide( pperp, density, m_q["Tperp"]);
    dg::blas1::pointwiseDivide( ppara, density, m_q["Tpara"]);

    //2. Solve for potential phi given density and temperature

    m_solvers.compute_phi( t, density, m_q["Tperp"], m_phi, m_p.penalize_wall,
        m_sources.get_wall(), m_p.penalize_sheath, m_parallel.get_sheath());

    m_solvers.compute_psi( t, density, m_q["Tperp"], m_phi,
        m_q["Psi0"], m_q["Psi1"], m_q["Psi2"], m_q["Psi3"]);

    timer.toc();
    accu += timer.diff();
    DG_RANK0 std::cout << "## Compute phi and psi               took "
                       << timer.diff()<<"s\t A: "<<accu<<"s\n";
    timer.tic( );
    // Given densities transform to fieldaligned grids:
    // (Guaranteed to compute "ST N")
    m_parallel.compute_staggered_densities( y, m_q);

    // Compute m_aparST if beta != 0
    if( m_p.beta != 0)
    {
        m_solvers.compute_aparST( t, m_q["ST N"], wST,
            m_aparST, true);
    }
    m_parallel.compute_apar( m_aparST, m_apar);
    m_old_apar.update( t, m_apar);
    timer.toc();
    accu += timer.diff();
    DG_RANK0 std::cout << "## Compute Apar and staggered N      took "
                       << timer.diff()<<"s\t A: "<<accu<<"s\n";
    timer.tic();
    for( unsigned s=0; s<m_p.num_species; s++)
        dg::blas1::axpby( 1., wST[s], -m_p.z[s]/m_p.mu[s], m_aparST, m_q["ST U"][s]);

    // Compute all the rest of parallel trafos
    m_parallel.compute_parallel_transformations( y, m_q);
    timer.toc();
    accu += timer.diff();
    DG_RANK0 std::cout << "## Compute Parallel transformations  took "
                       << timer.diff()<<"s\t A: "<<accu<<"s\n";
    timer.tic();
    // and the perp derivatives
    m_perp.update_derivatives( m_apar, m_dxapar, m_dyapar, y, m_q);
    m_perp.update_STderivatives( m_aparST, m_dxaparST, m_dyaparST, y, m_q);
    // Set perpendicular advection in yp
    // main species loop
    for( unsigned s=0; s<m_p.num_species; s++)
    {
        m_upToDate[s] = false;
        // Add perp dynamics
        m_perp.add_densities_advection(  s, m_apar, m_dxapar, m_dyapar,
            m_q, yp);
        m_perp.add_velocities_advection( s, m_aparST, m_dxaparST, m_dyaparST,
            m_q, yp);
        m_perp.add_densities_diffusion(  s, m_q, yp);
        m_perp.add_velocities_diffusion( s, m_q, yp);

        // Add parallel dynamics
        m_parallel.add_densities_advection(  s, y, m_q, yp);
        m_parallel.add_velocities_advection( s, y, m_q, yp);

        m_parallel.add_densities_diffusion( s, m_q, yp);
        m_parallel.add_velocities_diffusion( s, m_q, yp);

        m_parallel.add_sheath_neumann_terms( s, m_q, yp);
        m_parallel.add_sheath_velocity_terms( s, m_q, yp);

        // Add collisions
        m_collisions.add_coulomb_collisions( s, m_q, yp);
        m_collisions.add_lorentz_collisions( s, m_q, y, yp);

        // And sources
        m_sources.add_wall_terms( s, yp);
        m_sources.add_source_terms( s, m_perp, m_parallel, m_phi, m_q, y, yp );

        // Add penalization
        for( unsigned u=0; u<6; u++)
        {
            common::multiply_rhs_penalization( yp[u][s], m_p.penalize_wall,
                m_sources.get_wall(), m_p.penalize_sheath, m_parallel.get_sheath()); // F*(1-chi_w-chi_s)
            if( u == 3)
                // Apply to U, not W
                dg::blas1::pointwiseDot( -m_p.wall_rate, m_sources.get_wall(), m_q["ST U"][s],
                    -m_parallel.get_sheath_rate(), m_parallel.get_sheath(), m_q["ST U"][s], 1., yp[u][s]);
            else
                dg::blas1::pointwiseDot( -m_p.wall_rate, m_sources.get_wall(), y[u][s],
                    -m_parallel.get_sheath_rate(), m_parallel.get_sheath(), y[u][s], 1., yp[u][s]);
        }
    }

    timer.toc();
    accu += timer.diff();
    DG_RANK0 std::cout <<"## compute perp and para dynamics    took "
                       << timer.diff() << "s\t A: "<<accu<<"s\n";
}



} //namespace thermal
