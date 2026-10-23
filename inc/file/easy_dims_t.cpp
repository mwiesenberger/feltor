#include <iostream>
#include <filesystem>
#include <string>
#include <netcdf.h>
#include <cmath>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "catch2/catch_all.hpp"

#include "dg/algorithm.h"
#include "easy_dims.h"

// DEPRECATED
TEST_CASE( "Easy dims")
{
#ifdef WITH_MPI
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    //create a grid and some data
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {1, 1, 1});
#endif
    std::string filename = "dims.nc";
    INFO( "Write and check a couple of dimensions in "
                       << filename<<"\n");
    double x0 = 0., x1 = 2.*M_PI;
    dg::x::CartesianGrid3d grid( x0,x1,x0,x1,x0,x1,3,10,10,20
#ifdef WITH_MPI
    , comm
#endif
    );
    dg::ClonePtr<dg::x::aGeometry2d> grid2d_ptr = grid.perp_grid();

    //create NetCDF File
    int ncid = 0;
    dg::file::NC_Error_Handle err;
    DG_RANK0 err = nc_create( filename.data(), NC_NETCDF4|NC_CLOBBER, &ncid); //for netcdf4

    int dim_ids[4], tvarID;
    REQUIRE_NOTHROW(
        err = dg::file::define_dimensions( ncid, dim_ids, &tvarID, grid)
    );
    int dim2_ids[3];
    DG_RANK0// only master thread throws and knows
    {
        CHECK_THROWS_AS(
            err = dg::file::define_dimensions( ncid, dim2_ids, *grid2d_ptr),
            dg::file::NC_Error
        );
        bool exists = false;
        REQUIRE_NOTHROW(
            exists = dg::file::check_dimensions( ncid, dim2_ids, *grid2d_ptr)
        );
        INFO( "All dimensions exist?");
        CHECK( exists);
        REQUIRE_NOTHROW(
            err = dg::file::define_dimensions( ncid, dim2_ids, &tvarID, *grid2d_ptr, {"time", "x", "y"}, true)
        );
        INFO( "Check functions throw if existing dimensions do not match");
        int dimID, varID;
        dg::Grid1d g1d( 0., 1., 3, 5);
        CHECK( not dg::file::check_dimension( ncid, &dimID, g1d, "does-not-exist"));
        CHECK_THROWS_AS( dg::file::check_dimension( ncid, &dimID, g1d, "x"),
            std::runtime_error); // wrong length
        CHECK_THROWS_AS( dg::file::check_real_time<double>( ncid, "x", &dimID, &varID),
            std::runtime_error); // not unlimited
        CHECK_THROWS_AS( dg::file::check_real_time<float>( ncid, "time", &dimID, &varID),
            std::runtime_error); // wrong type
        CHECK( dg::file::check_real_time<double>( ncid, "time", &dimID, &varID));
        CHECK( dg::file::check_dimensions( ncid, dim2_ids, dg::Grid0d()));

        INFO( "Limited time");
        int lim_ids[3];
        dg::Grid2d g2d( 0., 1., 0., 1., 3, 4, 5);
        CHECK( dg::file::define_limtime_xy( ncid, lim_ids, 7, &varID, g2d,
            {"ltime", "ly", "lx"}) == NC_NOERR);
        size_t len;
        err = nc_inq_dimlen( ncid, lim_ids[0], &len);
        CHECK( len == 7);
        err = nc_inq_dimlen( ncid, lim_ids[1], &len);
        CHECK( len == 3*5);
        err = nc_inq_dimlen( ncid, lim_ids[2], &len);
        CHECK( len == 3*4);
    }

    DG_RANK0 err = nc_close(ncid);
    DG_RANK0 std::filesystem::remove( filename);
}
