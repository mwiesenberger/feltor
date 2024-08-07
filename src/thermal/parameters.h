#pragma once
#include <map>
#include <array>
#include <string>
#include <cmath>

#include "dg/file/json_utilities.h"

namespace thermal{
/// If you need more parameters, just go ahead and extend the list
struct Parameters
{
    unsigned n, Nx, Ny, Nz;

    unsigned itstp, maxout;
    double Tend, deltaT;
    unsigned cx, cy;

    std::vector<double> eps_pol;
    double jfactor;
    double eps_gamma, eps_ampere;
    unsigned stages;
    unsigned mx, my;
    double rk4eps;
    std::string interpolation_method;
    std::vector<double> nbc; // boundary condition for density

    double nu_ref, beta;

    unsigned diff_order;
    double nu_perp_n, nu_perp_u, nu_parallel_n;
    enum dg::direction diff_dir;
    std::string slope_limiter;

    double source_rate, nwall, uwall, wall_rate;

    enum dg::bc bcxN, bcyN, bcxU, bcyU, bcxP, bcyP, bcxA, bcyA;
    enum dg::direction pol_dir;
    std::string curvmode;
    std::string sheath_bc;
    std::string fci_bc;
    std::string output;
    bool symmetric, calibrate;
    bool penalize_wall, penalize_sheath;
    bool partitioned;
    //
    //
    unsigned num_species, num_trivial; // number of species and number of species with neglected mass ( electrons)
    double qlandau;
    std::vector<bool> neglect_mass;
    std::vector<int> z;
    std::vector<double> mu;
    std::vector<double> pi, kappa; // prefactors for Braginskii viscosity and conductivity
    std::vector<std::string> name; // name of species s

    Parameters() = default;
    Parameters( const dg::file::WrappedJsonValue& js) {
        //We need to check if a member is present
        n           = js["grid"].get("n", 3).asUInt();
        Nx          = js["grid"].get("Nx", 0).asUInt();
        Ny          = js["grid"].get("Ny", 0).asUInt();
        Nz          = js["grid"].get("Nz", 0).asUInt();
        partitioned = false;
        output      = js["output"].get( "type", "netcdf").asString();
        if( !("netcdf" == output) && !("glfw" == output))
            throw std::runtime_error( "Output type "+output+" not recognized!\n");
#ifdef WITHOUT_GLFW
        if( "glfw" == output)
            throw std::runtime_error( "Output type glfw not possible without glfw compiled!\n");
#endif
        cx = js["output"]["compression"].get(0u,1).asUInt();
        cy = js["output"]["compression"].get(1u,1).asUInt();

        stages      = js["elliptic"].get( "stages", 3).asUInt();
        eps_pol.resize(stages);
        eps_pol[0] = js["elliptic"]["eps_pol"].get( 0, 1e-6).asDouble();
        for( unsigned i=1;i<stages; i++)
        {
            eps_pol[i] = js["elliptic"][ "eps_pol"].get( i, 1).asDouble();
            eps_pol[i]*=eps_pol[0];
        }
        jfactor     = js["elliptic"].get( "jumpfactor", 1).asDouble();
        eps_gamma   = js["elliptic"].get( "eps_gamma", 1e-6).asDouble();
        eps_ampere  = js["elliptic"].get( "eps_ampere", 1e-6).asDouble();
        pol_dir = dg::str2direction(
                js["elliptic"].get("direction", "centered").asString() );


        mx          = js["FCI"]["refine"].get( 0u, 1).asUInt();
        my          = js["FCI"]["refine"].get( 1u, 1).asUInt();
        rk4eps      = js["FCI"].get( "rk4eps", 1e-6).asDouble();
        interpolation_method = js["FCI"].get("interpolation-method", "dg").asString();
        fci_bc      = js["FCI"].get( "bc", "along_field").asString();

        diff_order  = js["regularization"].get( "order", 2).asUInt();
        diff_dir    = dg::str2direction(
                js["regularization"].get( "direction", "centered").asString() );
        nu_perp_n   = js["regularization"].get( "nu_perp_n", 0.).asDouble();
        nu_perp_u   = js["regularization"].get( "nu_perp_u", 0.).asDouble();
        nu_parallel_n = js["regularization"].get( "nu_parallel_n", 0.).asDouble();
        slope_limiter = js["advection"].get("slope-limiter", "none").asString();
        if( (slope_limiter != "none") && (slope_limiter != "minmod")
             && (slope_limiter != "vanLeer")
                )
            throw std::runtime_error( "ERROR: advection : slope-limiter "+slope_limiter+" not recognized!\n");

        num_species = js["physical"]["species"].size();
        mu.resize( num_species);
        z.resize( num_species);
        neglect_mass.resize( num_species);
        num_trivial = 0;
        for( unsigned u=0; u<num_species; u++)
        {
            std::string name = js["physical"]["species"][u].asString();
            mu[u] = js["physical"]["mu"][u].asDouble();
            z[u] = js["physical"]["z"][u].asDouble();
            pi[u] = js["physical"]["pi"][u].asDouble();
            kappa[u] = js["physical"]["kappa"][u].asDouble();

            neglect_mass[u] = false;
            if( fabs(mu[u] ) < 1e-3)
            {
                neglect_mass[u] = true;
                num_trivial++;
            }
        }
        beta        = js["physical"].get( "beta", 0.).asDouble();
        nu_ref      = js["physical"].get( "collisionality", 0.).asDouble();
        qlandau     = js["physical"].get( "qlandau", 0.).asDouble();


        source_rate = 0.;
        if( js["source"].get("type", "zero").asString() != "zero")
            source_rate = js[ "source"].get( "rate", 0.).asDouble();
        sheath_bc = js["boundary"]["sheath"].get("type", "bohm").asString();
        if( (sheath_bc != "bohm") && (sheath_bc != "insulating") &&
                (sheath_bc != "none") && (sheath_bc != "wall"))
            throw std::runtime_error( "ERROR: Sheath bc "+sheath_bc+" not recognized!\n");

        bcxN = dg::str2bc(js["boundary"]["bc"][  "density"].get( 0, "").asString());
        bcyN = dg::str2bc(js["boundary"]["bc"][  "density"].get( 1, "").asString());
        nbc = 0.;
        if( bcxN == dg::DIR || bcxN == dg::DIR_NEU || bcxN == dg::NEU_DIR
            || bcyN == dg::DIR || bcyN == dg::DIR_NEU || bcyN == dg::NEU_DIR)
            nbc = js["boundary"]["bc"].get( "nbc", 1.0).asDouble();
            // assert that sum z[s] nbc[s] gives 0

        bcxU = dg::str2bc(js["boundary"]["bc"][ "velocity"].get( 0, "").asString());
        bcyU = dg::str2bc(js["boundary"]["bc"][ "velocity"].get( 1, "").asString());
        bcxP = dg::str2bc(js["boundary"]["bc"]["potential"].get( 0, "").asString());
        bcyP = dg::str2bc(js["boundary"]["bc"]["potential"].get( 1, "").asString());
        bcxA = dg::str2bc(js["boundary"]["bc"]["aparallel"].get( 0, "").asString());
        bcyA = dg::str2bc(js["boundary"]["bc"]["aparallel"].get( 1, "").asString());

        if( fci_bc == "along_field" || fci_bc == "perp")
        {
            if( bcxN != bcyN || bcxN == dg::DIR_NEU || bcxN == dg::NEU_DIR)
                throw std::runtime_error( "ERROR: density bc must be either dg::NEU or dg::DIR in both directions!\n");
            if( bcxU != bcyU || bcxU == dg::DIR_NEU || bcxU == dg::NEU_DIR)
                throw std::runtime_error( "ERROR: velocity bc must be either dg::NEU or dg::DIR in both directions!\n");
            if( bcxP != bcyP || bcxP == dg::DIR_NEU || bcxP == dg::NEU_DIR)
                throw std::runtime_error( "ERROR: potential bc must be either dg::NEU or dg::DIR in both directions!\n");
        }
        else if( fci_bc != "perp")
            throw std::runtime_error("Error! FCI bc '"+fci_bc+"' not recognized!\n");


        curvmode    = js["magnetic_field"].get( "curvmode", "toroidal").asString();
        penalize_wall = penalize_sheath = false;
        nwall = uwall = wall_rate = 0.;
        if( js["boundary"]["wall"].get("type","none").asString() != "none")
        {
            penalize_wall = js["boundary"]["wall"].get( "penalize-rhs",
                    false).asBool();
            wall_rate = js ["boundary"]["wall"].get( "penalization",
                    0.).asDouble();
            nwall = js["boundary"]["wall"].get( "nwall", 1.0).asDouble();
            uwall = js["boundary"]["wall"].get( "uwall", 0.0).asDouble();
        }
        if( sheath_bc != "none")
        {
            penalize_sheath = js["boundary"]["sheath"].get( "penalize-rhs",
                    false).asBool();
        }

        // Computing flags
        symmetric = calibrate = false;
        for( unsigned i=0; i<js["flags"].size(); i++)
        {
            std::string flag = js["flags"].get(i,"symmetric").asString();
            if( flag  == "symmetric")
                symmetric = true;
            else if( flag == "calibrate" )
            {
                if( output == "glfw")
                    throw std::runtime_error(
                            "Calibrate not possible with glfw output!\n");
                calibrate = true;
            }
            else
                throw std::runtime_error( "Flag "+flag+" not recognized!\n");
        }

        // output frequencies
        maxout = js["output"].get( "maxout", 0).asUInt();
        std::string output_mode = js["timestepper"].get(
                "output-mode", "Tend").asString();
        Tend = 0, deltaT = 0.;
        itstp       = js["output"].get("itstp", 0).asUInt();
        if( output_mode == "Tend")
        {
            Tend = js["timestepper"].get( "Tend", 1).asDouble();
            deltaT = Tend/(double)(maxout*itstp);
        }
        else if( output_mode == "deltaT")
        {
            deltaT = js["timestepper"].get( "deltaT", 1).asDouble()/(double)(itstp);
            Tend = deltaT*(double)(maxout*itstp);
        }
        else
        {
            throw std::runtime_error( "Error: Unrecognized timestepper: output-mode: '"
                               + output_mode);
        }
    }
};

}//namespace thermal
