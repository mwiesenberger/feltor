#ifndef _DG_PARAMETERS_ 
#define _DG_PARAMETERS_
#include <string>
#include "dg/enums.h"
#include "json/json.h"

namespace eule{
/**
 * @brief Provide a mapping between input file and named parameters
 */
struct Parameters
{
    unsigned n, Nx, Ny; 
    double dt; 
    unsigned n_out, Nx_out, Ny_out; 
    unsigned itstp, maxout;

    double eps_pol,  eps_gamma, eps_time;

    double mu[2];
    double tau[2];
    double mcv;
    double lx,ly;
    double nu_perp;
    
    double amp, sigma, posX, posY;

    double  nprofileamp, bgprofamp;
    unsigned init, iso, flrmode;
    enum dg::bc bc_x,bc_y; 

//     /**
//      * @brief constructor to make a const object
//      *
//      * @param v Vector from read_input function
//      */
//     Parameters( const std::vector< double>& v) {
//         n  = (unsigned)v[1]; 
//         Nx = (unsigned)v[2];
//         Ny = (unsigned)v[3];
//         dt = v[4];
//         n_out = v[5];
//         Nx_out = v[6];
//         Ny_out = v[7];
//         itstp = v[8];
//         maxout = v[9];
//         eps_pol = v[10];
//         eps_gamma = v[11];
//         eps_time = v[12];
//         mu[0] = -0.000272121;
//         mu[1] = 1.;
//         tau[0] = -1.;
//         tau[1]  = v[13];
//         mcv     = v[14];
//         nu_perp = v[15];
//         amp     = v[16];
//         sigma   = v[17];
//         posX    = v[18];
//         posY    = v[19];
//         nprofileamp = 0.;
//         bgprofamp   = 1.;
//         lx = v[20];
//         ly = v[21];
//         bc_x = map((int)v[22]);
//         bc_y = map((int)v[23]);
//         init = v[24];
//         iso =  v[25];
//         flrmode =  v[26];
//     }
    /**
     * @brief constructor to make a const object
     *
     * @param js json object
     */
    Parameters(const Json::Value& js) {
        n  = js["n"].asUInt();
        Nx = js["Nx"].asUInt();
        Ny = js["Ny"].asUInt();
        dt = js["dt"].asDouble();
        n_out  = js["n_out"].asUInt();
        Nx_out = js["Nx_out"].asUInt();
        Ny_out = js["Ny_out"].asUInt();
        itstp = js["itstp"].asUInt();
        maxout = js["maxout"].asUInt();
        
        eps_pol = js["eps_pol"].asDouble();
        eps_gamma = js["eps_gamma"].asDouble();
        eps_time = js["eps_time"].asDouble();
        mu[0]   = -0.000272121;
        mu[1]   = 1.;
        tau[0]  = -1.;
        tau[1]  = js["tau"].asDouble();
        mcv     = js["curvature"].asDouble();
        nu_perp = js["nu_perp"].asDouble();
        amp     = js["amplitude"].asDouble();
        sigma   = js["sigma"].asDouble();
        posX    = js["posX"].asDouble();
        posY    = js["posY"].asDouble();
        nprofileamp = 0.;
        bgprofamp   = 1.;
        lx = js["lx"].asDouble();
        ly =  js["ly"].asDouble();
        bc_x = dg::str2bc(js["bc_x"].asString());
        bc_y = dg::str2bc(js["bc_y"].asString());
        init = js["initmode"].asUInt();
        iso =  js["tempmode"].asUInt();
        flrmode =  js["flrmode"].asUInt();            
    }
    /**
     * @brief Display parameters
     *
     * @param os Output stream
     */
    void display( std::ostream& os = std::cout ) const
    {
        os << "Physical parameters are: \n"
            <<"     mu_e              = "<<mu[0]<<"\n"
            <<"     mu_i              = "<<mu[1]<<"\n"
            <<"     mcv               = "<<mcv<<"\n"
            <<"     El.-temperature:  = "<<tau[0]<<"\n"
            <<"     Ion-temperature:  = "<<tau[1]<<"\n"
            <<"     perp. Viscosity:  = "<<nu_perp<<"\n"
            <<"     eff grav./diss f. = "<<(1.+tau[1])*sigma*sigma*sigma*mcv*amp/(nu_perp*nu_perp)<<"\n"
            <<"     cst/dyn FLR (0/1) = "<<flrmode<<"\n"
            <<"     isothermal (0/1)  = "<<iso<<"\n";
        os  <<"Blob parameters are: \n"
            << "    amplitude:    "<<amp<<"\n"
            << "    width:        "<<sigma<<"\n"
            << "    posX:         "<<posX<<"\n"
            << "    posY:         "<<posY<<"\n";
        os << "Profile parameters are: \n"
            <<"     density profile amplitude:    "<<nprofileamp<<"\n"
            <<"     background profile amplitude: "<<bgprofamp<<"\n";
        os << "Algorithmic parameters are: \n"
            <<"     n  = "<<n<<"\n"
            <<"     Nx = "<<Nx<<"\n"
            <<"     Ny = "<<Ny<<"\n"
            <<"     dt = "<<dt<<"\n";
        os << "     Stopping for Polar CG:   "<<eps_pol<<"\n"
            <<"     Stopping for Gamma CG:   "<<eps_gamma<<"\n"
            <<"     Stopping for Time  CG:   "<<eps_time<<"\n";
        os << "Output parameters are: \n"
            <<"     n_out  =              "<<n_out<<"\n"
            <<"     Nx_out =              "<<Nx_out<<"\n"
            <<"     Ny_out =              "<<Ny_out<<"\n"
            <<"     Steps between output: "<<itstp<<"\n"
            <<"     Number of outputs:    "<<maxout<<"\n";
        os << "Box params: \n"
            <<"     lx  =              "<<lx<<"\n"
            <<"     ly  =              "<<ly<<"\n"; 
            displayBC( os, bc_x, bc_y);
        os << std::flush;//the endl is for the implicit flush 
    }
    private:
    dg::bc map( int i)
    {
        switch( i)
        {
            case(0): return dg::PER;
            case(1): return dg::DIR;
            case(2): return dg::DIR_NEU;
            case(3): return dg::NEU_DIR;
            case(4): return dg::NEU;
            default: return dg::PER;
        }
    }
    void displayBC( std::ostream& os, dg::bc bcx, dg::bc bcy) const
    {
        os << "Boundary conditions in x are: \n";
        switch( bcx)
        {
            case(0): os << "    PERIODIC";
                     break;
            case(1): os << "    DIRICHLET";
                     break;
            case(2): os << "    DIR_NEU";
                     break;
            case(3): os << "    NEU_DIR";
                     break;
            case(4): os << "    NEUMANN";
                     break;
        }
        os << "\nBoundary conditions in y are: \n";
        switch( bcy)
        {
            case(0): os << "    PERIODIC";
                     break;
            case(1): os << "    DIRICHLET";
                     break;
            case(2): os << "    DIR_NEU";
                     break;
            case(3): os << "    NEU_DIR";
                     break;
            case(4): os << "    NEUMANN";
                     break;
        }
        os <<"\n";
    }
};

}//namespace eule

#endif//_DG_PARAMETERS_

    

