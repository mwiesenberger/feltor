#pragma once


#include "dg/file/nc_utilities.h"
#include "parameters.h"
#include "thermal.h"

namespace thermal
{

using Feltor = thermal::Explicit< dg::x::CylindricalGrid3d, dg::x::IDMatrix,
        dg::x::DMatrix, dg::x::DVec>;

std::vector<dg::file::Record<void(dg::x::DVec&, Feltor&, unsigned i)>> restart3d_list = {
    {"restart_density", "gyro-centre density",
        []( dg::x::DVec& result, Feltor& f, unsigned i ) {
             dg::blas1::copy(f.restart_density(i), result);
        }
    },
    {"restart_velocity", "parallel velocity",
        []( dg::x::DVec& result, Feltor& f, unsigned i ) {
             dg::blas1::copy(f.restart_velocity(i), result);
        }
    },
    {"restart_pperp", "perpendicular pressure",
        []( dg::x::DVec& result, Feltor& f, unsigned i ) {
             dg::blas1::copy(f.restart_pperp(i), result);
        }
    },
    {"restart_ppara", "parallel pressure",
        []( dg::x::DVec& result, Feltor& f, unsigned i ) {
             dg::blas1::copy(f.restart_ppara(i), result);
        }
    },
    {"restart_aparallel", "parallel magnetic potential",
        []( dg::x::DVec& result, Feltor& f, unsigned i ) {
             dg::blas1::copy(f.restart_aparallel(), result);
        }
    }
};

std::array<std::vector<dg::x::DVec>,4> init_from_file( std::string file_name,
        const dg::x::CylindricalGrid3d& grid, const Parameters& p,
        double& time)
{
#ifdef WITH_MPI
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
    std::array<std::vector<dg::x::DVec>,4> y0;
    ///////////////////read in and show inputfile

    dg::file::NC_Error_Handle errIN;
    int ncidIN;
    errIN = nc_open( file_name.data(), NC_NOWRITE, &ncidIN);
    dg::file::WrappedJsonValue atts( dg::file::nc_attrs2json( ncidIN, NC_GLOBAL));
    dg::file::WrappedJsonValue jsIN = dg::file::string2Json(
        atts["inputfile"].asString(), dg::file::comments::are_forbidden);
    thermal::Parameters pIN( jsIN);
    DG_RANK0 std::cout << "# RESTART from file "<<file_name<< std::endl;
    DG_RANK0 std::cout << "#  file parameters:" << std::endl;
    DG_RANK0 std::cout << pIN.n<<" x "<<pIN.Nx<<" x "<<pIN.Ny<<" x "<<pIN.Nz
                <<" : symmetric "<<std::boolalpha<<pIN.symmetric<<std::endl;

    // Now read in last timestep
    dg::x::CylindricalGrid3d grid_IN( grid.x0(), grid.x1(), grid.y0(),
        grid.y1(), grid.z0(), grid.z1(),
        pIN.n, pIN.Nx, pIN.Ny, pIN.Nz, dg::DIR, dg::DIR, dg::PER
        #ifdef WITH_MPI
        , grid.communicator()
        #endif //WITH_MPI
        );
    // Theoretically we can change resolution
    dg::x::IHMatrix interpolateIN;
    dg::x::HVec transferIN;
    if( pIN.symmetric)
    {
        std::unique_ptr<dg::x::aGeometry2d> grid_perp( grid.perp_grid());
        interpolateIN = dg::create::interpolation( grid, *grid_perp);
        transferIN = dg::evaluate(dg::zero, *grid_perp);
    }
    else
    {
        interpolateIN = dg::create::interpolation( grid, grid_IN);
        transferIN = dg::evaluate(dg::zero, grid_IN);
    }

    dg::file::Reader<dg::x::CylindricalGrid3d> restart( ncidIN, grid, {"zr", "yr", "xr"});
    dg::file::Reader<dg::x::Grid0d> reader0d( ncidIN, {}, {"time"});
    /////////////////////Get time length and initial data///////////////////////////
    unsigned size_time = reader0d.size();
    reader0d.get( "time", time, size_time-1);
    DG_RANK0 std::cout << "# Current time = "<< time <<  std::endl;

    dg::x::HVec transferOUTvec = dg::evaluate( dg::zero, grid);
    dg::x::HVec apar = transferOUTvec;
    restart.get( "restart_aparallel", transferIN);
    dg::blas2::gemv( interpolateIN, transferIN, apar);
    for( unsigned s=0; s<pIN.num_species; s++)
    {
        for( unsigned k=0; k<restart3d_list.size()-1; k++) // skip aparallel
        {
            restart.get( restart3d_list[k].name + "_" + pIN.species[s].name, transferIN);
            dg::blas2::gemv( interpolateIN, transferIN, transferOUTvec);
            if( k == 1) // velocity: convert to W
                dg::blas1::axpby( 1., transferOUTvec, 1./p.species[s].mu,
                    apar, transferOUTvec);
            dg::assign( transferOUTvec, y0[k][i]);
        }
    }
    errIN = nc_close(ncidIN);
    /// ///////////////Now Construct initial fields ////////////////////////
    return y0;
}

}//namespace thermal
