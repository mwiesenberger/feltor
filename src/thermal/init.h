#pragma once

#include "dg/file/json_utilities.h"
#include "init_from_file.h"
#include "../feltor/init.h"

namespace thermal
{

/* The purpose of this file is to provide an interface for custom initial conditions and
 * source profiles.  Just add your own to the relevant map below.
 * @note y0[3,4,5] havec to be staggered velocity, qperp and qpara
 */

std::array<std::vector<dg::x::DVec>,6> initial_conditions(
    Explicit<dg::x::CylindricalGrid3d, dg::x::IDMatrix, dg::x::DMatrix, dg::x::DVec>& thermal,
    const dg::x::CylindricalGrid3d& grid, const thermal::Parameters& p,
    const dg::geo::TokamakMagneticField& mag,
    const dg::geo::TokamakMagneticField& unmod_mag,
    dg::file::WrappedJsonValue jsspecies, // input["species"]
    double & time, dg::geo::CylindricalFunctor& sheath_coordinate )
{
#ifdef WITH_MPI
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
    // First, let's scan for restart
    for( unsigned s=0; s<p.num_species; s++)
    {
        dg::file::WrappedJsonValue js = jsspecies[s]["init"];
        std::string type = js["type"].asString();
        if( "fields" == type)
            continue;
        else if( "restart" == type)
        {
            std::string file = js["file"].asString();
            return init_from_file( file, grid, p, time);
        }
        else
            throw dg::Error(dg::Message()<< "Invalid initial condition "<<type<<"\n");
    }
    time = 0;
    std::array<std::vector<dg::x::DVec>,6> y0;
    for( unsigned i=0; i<6; i++)
    for( unsigned s=0; s<p.num_species; s++)
        y0[i][s] = dg::construct<dg::x::DVec>( dg::evaluate( dg::zero, grid));
    unsigned num_quasineutral_n = 0, idx_quasineutral_n;
    // Init physical n, pperp, ppara
    for( unsigned s=0; s<p.num_species; s++)
    {
        dg::file::WrappedJsonValue js = jsspecies[s]["init"];
        DG_RANK0  std::cout << "# Initializing species "<<p.name[s]<<"\n";

        const std::array<std::string,3> eqs = { "density", "pperp", "ppara"};
        // Think of this as initialising the **physical** quantities
        for( std::string field : eqs)
        {
            std::string ntype = js[field].get("type", "zero").asString();
            DG_RANK0 std::cout << "# Initialize "<<field<<" with "<<ntype << std::endl;
            if ( "const" == ntype)
            {
                double nbg = js[field].get("background", 0.1).asDouble();
                y0[eq][s] = dg::construct<dg::x::DVec>(
                        dg::evaluate( dg::CONSTANT( nbg), grid));
            }
            else if( "ne" == ntype )
            {
                double nbg = 0.;
                dg::x::HVec ntilde  = detail::make_ntilde(  thermal, grid, mag,
                        js[field]["ntilde"]);
                dg::x::HVec profile = detail::make_profile( grid, mag,
                        js[field]["profile"], nbg);
                dg::x::HVec damping = detail::make_damping( grid, unmod_mag,
                        js[field]["damping"]);
                dg::x::HVec density = profile;
                dg::blas1::subroutine( [nbg]( double profile, double
                    ntilde, double damping, double& density)
                        { density = (profile+ntilde-nbg)*damping+nbg;},
                        profile, ntilde, damping, density);
                dg::assign( density, y0[eq][s]);
            }
            else if ( "density" == field and "quasineutral" == ntype)
            {
                num_quasineutral_n ++ ;
                idx_quasineutral_n = s;
            }
            else
                throw dg::Error(dg::Message()<< "Invalid "<<field<<" initial condition "<<ntype<<"\n");
        }
    }
    // Check quasineutrality
    if( num_quasineutral_n > 1)
        throw std::runtime_error( "There cannot be more than one quasineutral density");
    if( num_quasineutral_n == 0)
        ;//TODO Check that species indeed add up to 0
    else
    {
        for( unsigned s=0; s<p.num_species and s != idx_quasineutral_n; s++)
            dg::blas1::axpby( -p.z[s]/ p.z[idx_quasineutral_n], y0[0][s], 1.,
                y0[0][idx_quasineutral_n]);
    }
    // Now transform to gyro-centre density
    for( unsigned s=0; s<p.num_species; s++)
    {
        dg::file::WrappedJsonValue js = jsspecies[s]["init"];
        std::string invert = js["density"].get("invert", "none").asString();
        if( "none" == invert)
            continue; // N_s = n_s
        //TODO Fix lwl
        //else if( "lwl" == invert)
        //{
        //    //add FLR correction -0.5*tau*mu*Delta n_e
        //    dg::blas1::transform(src, m_temp0, dg::PLUS<double>(-m_p.nbc));
        //    dg::blas2::symv( 0.5*m_p.tau[1]*m_p.mu[1],
        //        m_lapperpN, m_temp0, 1.0, target);
        //    double minimalni = dg::blas1::reduce( y0[0][1], 1e10,
        //            thrust::minimum<double>());
        //    DG_RANK0 std::cerr << "# Minimum Ni value "<<minimalni<<std::endl;
        //    //actually we should always invert according to Markus
        //    //because the dG direct application is supraconvergent
        //    std::string ptype = js["potential"].get("type", "zero_pol").asString();
        //    DG_RANK0 std::cout << "# Initialize potential with "<<ptype << std::endl;
        //    if( minimalni <= 0.0)
        //    {
        //        throw dg::Error(dg::Message()<< "ERROR: invalid initial condition. Increase value for alpha since now the ion gyrocentre density is negative!\n"
        //            << "Minimum Ni value "<<minimalni);
        //    }
        //}
        else
            throw dg::Error(dg::Message()<< "Invalid density invert type "<<invert<<"\n");

    }

    unsigned num_quasineutral_u = 0, idx_quasineutral_u;
    // Init physical u, qperp, qpara
    for( unsigned s=0; s<p.num_species; s++)
    {
        // init (staggered) velocity and thus canonical W
        std::string utype = js["velocity"].get("type", "zero").asString();
        DG_RANK0 std::cout << "# Initialize velocity with "<<utype << std::endl;
        if( "zero" == utype)
            ; // velocity already is zero
        else if( "linear_cs" == utype )
        {
            // TODO fix sheath init, assume toroidal alignment ot T in sheath
            //std::string uprofile = js["velocity"].get("profile", "linear_cs").asString();

            //dg::blas1::evaluate( m_temp, dg::equals(),
            //    [ zs = m_p.z[s], ms = m_p.mu[s]] DG_DEVICE( double ts, double te)
            //    {
            //        return sqrt( ( ts + zs * te )/ms);
            //    }, tparaST[k], tparaST[0])
            //std::unique_ptr<dg::x::aGeometry2d> perp_grid_ptr( grid.perp_grid());
            //dg::x::HVec coord2d = dg::pullback( sheath_coordinate,
            //        *perp_grid_ptr);
            //dg::x::HVec ui;
            //dg::assign3dfrom2d( coord2d, ui, grid);
            //dg::blas1::scal( ui, sqrt( 1.0+p.tau[1]));
            //dg::assign( ui, y0[1][1]); // Wi = Ui
        }
        else if( "quasineutral" == utype)
        {
            num_quasineutral_u ++;
            idx_quasineutral_u = s;
        }
        else
            throw dg::Error(dg::Message()<< "Invalid velocity initial condition "<<utype<<"\n");

        std::string qperptype = js["qperp"].get("type", "zero").asString();
        DG_RANK0 std::cout << "# Initialize qperp with "<<utype << std::endl;
        if( "zero" == qperptype)
            ; // qperp already is zero
        else
            throw dg::Error(dg::Message()<< "Invalid qperp initial condition "<<qperptype<<"\n");

        std::string qparatype = js["qpara"].get("type", "zero").asString();
        DG_RANK0 std::cout << "# Initialize qperp with "<<utype << std::endl;
        if( "zero" == qparatype)
            ; // qpara
        else
            throw dg::Error(dg::Message()<< "Invalid qpara initial condition "<<qparatype<<"\n");
    }
    // Init quasineutral species
    if( num_quasineutral_u > 1)
        throw std::runtime_error( "There cannot be more than one quasineutral velocity");
    if( num_quasineutral_u == 0)
        ;//TODO Check that species indeed add up to 0
    else
    {
        // Assume densities are toroidally aligned (so nST = n)
        for( unsigned s=0; s<p.num_species and s != idx_quasineutral_u; s++)
            dg::blas1::evaluate( y0[3][idx_quasineutral_u], dg::plus_equals(),
                [ zk = p.z[idx_quasineutral_u], zs = p.z[s]] DG_DEVICE(
                    double nk, double ns, double us)
                {
                    return -zs * ns * us / (zk * nk);
                }, y0[0][idx_quasineutral_u], y0[0][s], y0[3][s])
    }
    // Since total current is 0, we have A_\parallel = 0

    return y0;
};

std::array<std::vector<dg::x::HVec>,3> source_profiles(
    Explicit<dg::x::CylindricalGrid3d, dg::x::IDMatrix, dg::x::DMatrix, dg::x::DVec>& thermal,
    bool& fixed_profile,        //indicate whether a profile should be forced (yes or no)
    std::array<std::vector<dg::x::HVec>,3>& profiles,    // if fixed_profile is yes you need to construct something here, if no then you can ignore the parameter; if you construct something it will show in the output file
    const dg::x::CylindricalGrid3d& grid,
    const dg::geo::TokamakMagneticField& mag,
    const dg::geo::TokamakMagneticField& unmod_mag,
    dg::file::WrappedJsonValue jsspecies, // input["species"]
    std::vector<double>& minne,
    std::vector<double>& minrate,
    std::vector<double>& minalpha
    )
{
    unsigned num_species = jsspecies.size();
    minne.resize( num_species);
    minrate.resize( num_species);
    minalpha.resize( num_species);
    std::array<std::vector<dg::x::HVec>,3> source;
    for( unsigned s=0; s<num_species; s++)
    {
        dg::file::WrappedJsonValue js = jsspecies[s]["source"];
        minne[s] = js["density"].get("minne", 0.).asDouble();
        minrate[s] =  minalpha[s] = 0;
        if( minne != 0)
        {
            minrate[s] = js["density"].get("minrate", 1.).asDouble();
            minalpha[s] = js["density"].get("minalpha", 0.05).asDouble();
        }
        if( minrate[s] != minrate[0])
            throw dg::Error(dg::Message()<< "Minrate must be equal among species "<<minrate[s]<<" vs "<<minrate[0]<<"!\n");

        const std::array<std::string,3> eqs = { "density", "pperp", "ppara"};
        // Think of this as initialising the **physical** quantities
        for( unsigned u=0; u<3; u++)
        {
            source[u][s] = dg::evaluate( dg::zero, grid);

            std::string type  = js[eqs[u]].get( "type", "zero").asString();

            profiles[u][s] = dg::evaluate( dg::zero, grid);
            if( "zero" == type)
            {
                ;
            }
            else if( "fixed_profile" == type)
            {
                fixed_profile = true;
                double nbg = 0;
                ne_profile = detail::make_profile(grid, mag, js["profile"], nbg);
                source = detail::make_damping( grid, unmod_mag, js["damping"]);
            }
            else if("influx" == type)
            {
                fixed_profile = false;
                double nbg = 0.;
                source  = detail::make_ntilde( thermal, grid, mag, js["ntilde"]);
                ne_profile = detail::make_profile( grid, mag, js["profile"], nbg);
                dg::x::HVec damping = detail::make_damping( grid, unmod_mag, js["damping"]);
                dg::blas1::subroutine( [nbg]( double& profile, double& ntilde, double
                            damping) {
                            ntilde  = (profile+ntilde-nbg)*damping+nbg;
                            profile = (profile-nbg)*damping +nbg;
                        },
                        ne_profile, source, damping);
            }
            else if( "torpex" == type)
            {
                fixed_profile = false;
                double rhosinm = 0.98 / mag.R0();
                double rhosinm2 = rhosinm*rhosinm;
                source = dg::pullback( detail::TorpexSource(
                0.98/rhosinm, -0.02/rhosinm, 0.0335/rhosinm, 0.05/rhosinm, 565*rhosinm2),
                    grid);
            }
            else if( "tcv" == type)
            {
                fixed_profile = false;
                const double R_0 = 1075., Z_0 = -10.;
                const double psip0 = mag.psip()( R_0, Z_0);
                const double sigma = 9.3e-3*psip0/0.4;
                source = dg::pullback(
                    dg::compose( dg::GaussianX( psip0, sigma, 1.),  mag.psip() ), grid);
                dg::blas1::pointwiseDot( detail::xpoint_damping(grid,mag),
                       source, source);
            }
            else
                throw dg::Error(dg::Message()<< "Invalid source type "<<type<<"\n");
        }
    }
    return source;
};

} //namespace thermal
