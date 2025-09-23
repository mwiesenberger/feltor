#pragma once

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "perpendicular.h"
#include "parallel.h"
#include "parameters.h"

namespace thermal
{
// Sources and Damping
template< class Geometry, class IMatrix, class Matrix, class Container >
struct Sources
{
    Sources( const Geometry&, thermal::Parameters,
        dg::geo::TokamakMagneticField, dg::file::WrappedJsonValue);
    void add_wall_terms( unsigned s,  std::array<std::vector<Container>,6>& yp);
    void add_source_terms(
        unsigned s,
        PerpDynamics<Geometry, IMatrix, Matrix, Container>& perp,
        ParallelDynamics<Geometry, IMatrix, Matrix, Container>& parallel,
        const Container& phi,
        const std::map<std::string, std::vector<Container>>& q,
        const std::array<std::vector<Container>,6>& y,
        std::array<std::vector<Container>,6>& yp);
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
            m_profile[u].resize( m_p.num_species);
            m_source[u].resize( m_p.num_species);
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
        dg::assign( wall, m_wall);
    }
    const Container& get_wall() const{
        return m_wall;
    }
    const Container& get_source(unsigned u, unsigned s) const{
        return m_source[u][s];
    }
    const Container& get_source_prof(unsigned u, unsigned s) const{
        return m_profile[u][s];
    }
    private:
    const thermal::Parameters m_p;
    Container m_temp0, m_temp1;
    Container m_wall;
    std::array<std::vector<Container>,3> m_ss; // source terms for all species
    std::array<std::vector<Container>,3> m_source, m_profile; // (physical) source terms for all species

    std::array<std::vector<double>,3> m_source_rate;
    bool m_fixed_profile = false;

    std::vector<double> m_minne, m_minalpha;
    double m_minrate  = 0.;

};
template<class Grid, class IMatrix, class Matrix, class Container>
Sources<Grid, IMatrix, Matrix, Container>::Sources( const Grid& g,
    thermal::Parameters p, dg::geo::TokamakMagneticField,
    dg::file::WrappedJsonValue
    ): m_p(p)
{
    dg::assign( dg::evaluate( dg::zero, g), m_temp0 );
    m_temp1 = m_temp0;
    for( int i=0; i<3; i++)
        m_ss[i].resize( m_p.num_species, m_temp0);
    m_source = m_profile = m_ss;
}


template<class Geometry, class IMatrix, class Matrix, class Container>
void Sources<Geometry, IMatrix, Matrix, Container>::add_source_terms(
    unsigned s,
    PerpDynamics<Geometry, IMatrix, Matrix, Container>& perp,
    ParallelDynamics<Geometry, IMatrix, Matrix, Container>& parallel,
    const Container& phi,
    const std::map<std::string, std::vector<Container>>& q,
    const std::array<std::vector<Container>,6>& y,
    std::array<std::vector<Container>,6>& yp)
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
        perp.transform_density_pperp( m_p.mu[s], m_p.z[s], m_profile[0][s],
            m_profile[1][s], phi, m_ss[0][s], m_ss[1][s]);
        dg::blas1::copy( m_profile[2][s], m_ss[2][s]);
        for( unsigned u=0; u<3; u++)
        {
            if( m_source_rate[u][s] != 0.0)
            {
                dg::blas1::subroutine(
                    [] DG_DEVICE ( double& result, double ne, double profne,
                        double source, double source_rate){
                        result = source_rate*source*(profne - ne);
                        },
                    m_ss[u][s], y[u][s], m_ss[u][s], m_source[u][s], m_source_rate[u][s]);
            }
            else
                dg::blas1::copy( 0., m_ss[u][s]);
        }
    }
    else // influx
    {
        for( unsigned u=0; u<3; u++)
        {
            if( m_source_rate[u][s] != 0.0)
            {
                dg::blas1::axpby( m_source_rate[u][s], m_source[u][s], 0., m_ss[u][s]);
            }
            else
                dg::blas1::copy( 0., m_ss[u][s]);
        }
        // if influx transform sources last
        perp.transform_density_pperp( m_p.mu[s], m_p.z[s], m_ss[0][s],  m_ss[1][s], phi,
            m_ss[0][s], m_ss[1][s]);
    }
    //Add all to the right hand side
    for( unsigned u=0; u<3; u++)
        dg::blas1::axpby( 1., m_ss[u][s], 1.0, yp[u][s]);

    // TODO maybe this can be avoided by re-computing the density source on staggered grid...
    parallel.add_velocity_source_term( s, m_ss[0], q["ST N"], q["ST U"], yp);
}

template<class Geometry, class IMatrix, class Matrix, class Container>
void Sources<Geometry, IMatrix, Matrix, Container>::add_wall_terms(
    unsigned s,
    std::array<std::vector<Container>,6>& yp)
{
    // add wall boundary conditions
    if( m_p.wall_rate != 0)
    {
        dg::blas1::axpby( +m_p.wall_rate*m_p.nwall[s], m_wall, 1., yp[0][s] );
        dg::blas1::axpby( +m_p.wall_rate*m_p.twall*m_p.nwall[s], m_wall, 1., yp[1][s] );
        dg::blas1::axpby( +m_p.wall_rate*m_p.twall*m_p.nwall[s], m_wall, 1., yp[2][s] );
        dg::blas1::axpby( +m_p.wall_rate*m_p.uwall, m_wall, 1., yp[3][s] );
        dg::blas1::axpby( +m_p.wall_rate*m_p.qwall, m_wall, 1., yp[4][s] );
        dg::blas1::axpby( +m_p.wall_rate*m_p.qwall, m_wall, 1., yp[5][s] );
    }
}

} //namespace thermal
