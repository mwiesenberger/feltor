#include <iostream>
#include <memory>
#include <mpi.h>

#include "dg/algorithm.h"
#include "mpi_curvilinear.h"
#include "polar.h"
#include "catch2/catch_all.hpp"

TEST_CASE( "MPI curvilinear grids")
{
    const dg::geo::PolarGenerator generator( 1, 2);
    SECTION( "Global 3d geometry")
    {
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0,0}, {false,true,true});
        // The z boundaries differ from the boundaries in x
        const dg::RealMPIGrid1d<double> gz( 0.5, 1.5, 1, 4, dg::PER,
                dg::mpi_cart_split_as<3>(comm)[2]);
        const dg::geo::CurvilinearProductMPIGrid3d grid( generator,
                {3,4,dg::DIR}, {3,8,dg::PER}, gz, comm);
        const std::unique_ptr<dg::aGeometry3d> global( grid.global_geometry());
        for( unsigned u=0; u<3; u++)
        {
            INFO( "Direction "<<u);
            CHECK( global->get_p()[u] == grid.global().get_p()[u]);
            CHECK( global->get_q()[u] == grid.global().get_q()[u]);
            CHECK( global->get_n()[u] == grid.global().get_n()[u]);
            CHECK( global->get_N()[u] == grid.global().get_N()[u]);
            CHECK( global->get_bc()[u] == grid.global().get_bc()[u]);
        }
        CHECK( global->z0() == 0.5);
        CHECK( global->z1() == 1.5);
    }
    SECTION( "Global 2d geometry")
    {
        MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {false,true});
        const dg::geo::CurvilinearMPIGrid2d grid( generator, 3, 4, 8, dg::DIR, dg::PER, comm);
        const std::unique_ptr<dg::aGeometry2d> global( grid.global_geometry());
        for( unsigned u=0; u<2; u++)
        {
            INFO( "Direction "<<u);
            CHECK( global->get_p()[u] == grid.global().get_p()[u]);
            CHECK( global->get_q()[u] == grid.global().get_q()[u]);
            CHECK( global->get_n()[u] == grid.global().get_n()[u]);
            CHECK( global->get_N()[u] == grid.global().get_N()[u]);
            CHECK( global->get_bc()[u] == grid.global().get_bc()[u]);
        }
    }
}
