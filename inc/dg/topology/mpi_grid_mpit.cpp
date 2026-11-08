#include <iostream>
#include <sstream>
#include <mpi.h>
#include "dg/backend/mpi_init.h"
#include "mpi_grid.h"
#include "mpi_evaluation.h"

#include "catch2/catch_all.hpp"

static int function( double x, double y)
{
    return int(10000*x*y);
}
// TODO Test what happens if # processes is larger than N?
TEST_CASE( "MPI Grid")
{
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    unsigned nx = 2, ny = 3,  Nx = 8, Ny = 10;
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {1,1});
    MPI_Comm comm1d = dg::mpi_cart_create( MPI_COMM_WORLD, {0}, {1});
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // Assume that evaluate works
    dg::MPIGrid2d g2d( {0.,0.}, {1.,1.}, {nx,ny}, {Nx,Ny}, {dg::PER, dg::PER},
            dg::mpi_cart_split_as<2>(comm));
    auto local_vec = dg::evaluate( function, g2d);
    auto global_vec = dg::evaluate( function, g2d.global());
    SECTION( "global2local")
    {
        auto local_vec2 = dg::global2local( global_vec, g2d);
        int local_size  = g2d.local_size();
        for( int i=0; i<local_size; i++)
            CHECK( local_vec.data()[i] == local_vec2.data()[i]);
    }

    SECTION("global2localIdx")
    {
        int global_size = g2d.size();
        for ( int i=0; i<global_size; i++)
        {
            int lIdx = 0, pid = 0;
            g2d.global2localIdx( i, lIdx, pid);
            if( pid == rank)
            {
                CHECK( global_vec[i] == local_vec.data()[lIdx]);
                CHECK( lIdx < (int)local_vec.data().size());
            }
        }
    }

    SECTION( "Processes agree on indices")
    {
        int global_size = g2d.size();
        for ( int i=0; i<global_size; i++)
        {
            int lIdx = 0, pid = 0;
            g2d.global2localIdx( i, lIdx, pid);
            std::vector<int> locals(size), localp(size);
            // Send everything to root
            MPI_Gather( &lIdx, 1, MPI_INT, &locals[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(  &pid, 1, MPI_INT, &localp[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(rank==0)
            {
                for( int k=0; k<size; k++)
                    CHECK( locals[k] == lIdx);
                for( int k=0; k<size; k++)
                    CHECK( localp[k] == pid);
            }
            int gIdx = 0;
            g2d.local2globalIdx( lIdx, pid, gIdx);
            INFO("Rank "<<rank<<" pid "<<pid <<" lIdx "<<lIdx<< " gIdx "<<gIdx<<" "<<i);
            CHECK( gIdx == i);
        }
    }
    SECTION( "Shape is global")
    {
        std::array<unsigned,2> shape = { nx*Nx, ny*Ny};
        CHECK( g2d.shape(0) == shape[0]);
        CHECK( g2d.shape(1) == shape[1]);
        CHECK( g2d.get_shape() == shape);
    }
    SECTION( "Global grids larger than the int range")
    {
        const unsigned big = 65536;
        dg::MPIGrid2d gbig( {0.,0.}, {1.,1.}, {1u,1u}, {big,big},
            {dg::PER, dg::PER}, dg::mpi_cart_split_as<2>(comm));
        CHECK( gbig.size() == size_t(1) << 32);
        const int int_max = std::numeric_limits<int>::max();
        int lIdx = -1, pid = -1, gIdx = -1;
        // All int indices are valid global indices
        CHECK( gbig.global2localIdx( int_max, lIdx, pid));
        // The last point has a global index that does not fit into int
        int cart_rank, cart_size;
        MPI_Comm_rank( gbig.communicator(), &cart_rank);
        MPI_Comm_size( gbig.communicator(), &cart_size);
        size_t last = gbig.local_size()-1;
        if( cart_rank == cart_size-1 and last <= size_t(int_max))
            CHECK( not gbig.local2globalIdx( int(last), cart_rank, gIdx));
        // Small global indices keep working
        if( cart_rank == 0)
        {
            CHECK( gbig.local2globalIdx( 0, cart_rank, gIdx));
            CHECK( gIdx == 0);
        }
    }
    SECTION( "local2globalIdx rejects indices outside the local vector")
    {
        for( int pid=0; pid<size; pid++)
        {
            // Every process can compute the local size of pid
            int lsize = 0, gIdx = 0;
            for( int i=0; i<(int)g2d.size(); i++)
            {
                int lIdx = 0, p = 0;
                g2d.global2localIdx( i, lIdx, p);
                if( p == pid)
                    lsize++;
            }
            INFO( "Rank "<<rank<<" pid "<<pid<<" local size "<<lsize);
            CHECK( g2d.local2globalIdx( lsize-1, pid, gIdx));
            CHECK( not g2d.local2globalIdx( lsize, pid, gIdx));
            CHECK( not g2d.local2globalIdx( -1, pid, gIdx));
        }
    }
    SECTION( "Display")
    {
        std::stringstream ss;
        g2d.display( ss);
        // global and local grid parameters
        CHECK( ss.str().find( "GLOBAL GRID") != std::string::npos);
        CHECK( ss.str().find( "    N  = "+std::to_string(Nx)) != std::string::npos);
        CHECK( ss.str().find( "    N  = "+std::to_string(g2d.local().N(0)))
                != std::string::npos);
        CHECK( ss.str().find( "LOCAL GRID") != std::string::npos);
        CHECK( ss.str().find( "LOCAL GRID") < ss.str().rfind( "    N  = "));
    }
    SECTION( "Start and count")
    {
        int gIdx;
        g2d.local2globalIdx( 0, rank, gIdx);
        auto start = g2d.start();
        int sIdx = start[0]*g2d.shape(0) + start[1];
        INFO( "Rank "<<rank<<" sIdx "<<sIdx<<" gIdx "<<gIdx);
        CHECK( sIdx == gIdx);
    }

    SECTION( "Local and global boundaries")
    {
        // This reproduces a bug in a TCV feltor3d run
        double boxscaleZm = 3.8, boxscaleZp = 1.85;
        const double scale = 1000.5340440422314;
        double R0 = 0.9*scale;
        double inverseaspectratio = 0.2777777777777778;
        double a = R0*inverseaspectratio;
        const double Zmin = -boxscaleZm*a;
        const double Zmax = +boxscaleZp*a;
        dg::MPIGrid1d g1d( Zmin, Zmax, 3, 80, comm1d );
        if( rank == 0)
            CHECK( g1d.local().p() == g1d.p());
        if( rank == size-1)
            CHECK( g1d.local().q() == g1d.q());
        // The setters must give the same local grid as the constructor
        auto check_local = [&]( const dg::MPIGrid1d& g)
        {
            dg::MPIGrid1d g_new( g.p(), g.q(), g.n(), g.N(), comm1d);
            INFO( "Rank "<<rank<<" n "<<g.n()<<" N "<<g.N());
            CHECK( g.local().p() == g_new.local().p());
            CHECK( g.local().q() == g_new.local().q());
            CHECK( g.local().N() == g_new.local().N());
            CHECK( g.local().n() == g_new.local().n());
        };
        g1d.set_pq( {0.}, {1.});
        check_local( g1d);
        g1d.set_pq( {Zmin}, {Zmax});
        check_local( g1d);
        g1d.set( 2, 81);
        check_local( g1d);
        g1d.set( 3, 79);
        check_local( g1d);
    }
}
