#include <iostream>

#include "../blas.h"
#include "transform.h"
#include "refined_gridX.h"
#include "dlt.h"

#include "catch2/catch_all.hpp"

// Cell boundaries reconstructed from the abscissas of one line of a grid
static std::vector<double> cell_boundaries( const std::vector<double>& abs,
    unsigned n, double x0)
{
    std::vector<double> xi = dg::DLT<double>::abscissas(n);
    std::vector<double> b( 1, x0);
    for( unsigned i=0; i<abs.size()/n; i++)
    {
        double width = 2.*(abs[i*n+n-1]-abs[i*n])/(xi[n-1]-xi[0]);
        b.push_back( (abs[i*n+n-1]+abs[i*n])/2. + width/2.);
    }
    return b;
}

// Check that a refined X grid covers the same domain as the unrefined
// grid with the X-point on a cell boundary; return the cell boundaries
static std::array<std::vector<double>,2> check_refinedX(
    const dg::aRefinementX2d& ref, double fx, double fy, dg::bc bcx)
{
    const double x0 = 0, x1 = 1, y0 = 0, y1 = 2;
    const unsigned n = 3, Nx = 8, Ny = 8;
    dg::GridX2d g( x0, x1, y0, y1, fx, fy, n, Nx, Ny, bcx, dg::NEU);
    dg::CartesianRefinedGridX2d g_f( ref, x0, x1, y0, y1, fx, fy, n, Nx,
            Ny, bcx, dg::NEU);
    std::vector<dg::HVec> map = g_f.map();
    REQUIRE( map[0].size() == g_f.size());
    REQUIRE( map[1].size() == g_f.size());
    unsigned NX = n*g_f.Nx();
    std::vector<double> ax( map[0].begin(), map[0].begin()+NX), ay;
    for( unsigned i=0; i<n*g_f.Ny(); i++)
        ay.push_back( map[1][i*NX]);
    std::vector<double> bx = cell_boundaries( ax, n, x0);
    std::vector<double> by = cell_boundaries( ay, n, y0);
    CHECK( fabs( bx.back() - x1) < 1e-13);
    CHECK( fabs( by.back() - y1) < 1e-13);
    // The X-point lies on the same cell boundaries
    CHECK( fabs( bx[g_f.inner_Nx()] - (x0 + g.inner_Nx()*g.hx())) < 1e-13);
    CHECK( fabs( by[g_f.outer_Ny()] - (y0 + g.outer_Ny()*g.hy())) < 1e-13);
    CHECK( fabs( by[g_f.outer_Ny()+g_f.inner_Ny()]
                - (y1 - g.outer_Ny()*g.hy())) < 1e-13);
    double area = dg::blas1::dot( dg::create::volume( g_f), 1.);
    CHECK( fabs( area - (x1-x0)*(y1-y0)) < 1e-13);
    return {bx, by};
}

// Number of cells of given width in a list of cell boundaries
static unsigned num_cells( const std::vector<double>& b, double width)
{
    unsigned num = 0;
    for( unsigned i=0; i+1<b.size(); i++)
        if( fabs( b[i+1]-b[i] - width) < 1e-13)
            num++;
    return num;
}

TEST_CASE( "Refinement X")
{
    SECTION( "Identity X Refinement")
    {
        // The identity refinement reproduces the X-point grid
        dg::IdentityXRefinement id;
        dg::CartesianRefinedGridX2d g2d_f( id, 0., 1., 0., 2., 0.25, 0.25,
                3, 8, 8, dg::DIR, dg::NEU);
        dg::CartesianGridX2d g2d( 0., 1., 0., 2., 0.25, 0.25, 3, 8, 8,
                dg::DIR, dg::NEU);
        REQUIRE( g2d_f.size() == g2d.size());
        std::vector<dg::HVec> map = g2d_f.map();
        std::vector<dg::HVec> sol = g2d.map();
        dg::HVec vol_f = dg::create::volume( g2d_f);
        dg::HVec vol = dg::create::volume( g2d);
        for( unsigned i=0; i<g2d.size(); i++)
        {
            INFO( "Index "<<i);
            CHECK( map[0][i] == sol[0][i]);
            CHECK( map[1][i] == sol[1][i]);
            CHECK( vol_f[i] == vol[i]);
        }
    }
    SECTION( "Equidist X Refinement")
    {
        // Each refined cell is divided into add+1 equidistant cells
        for( unsigned howm : {1u,2u})
        {
            INFO( "howmany "<<howm);
            dg::EquidistXRefinement ref( 2, 1, howm, howm);
            auto b = check_refinedX( ref, 0.25, 0.25, dg::DIR);
            CHECK( num_cells( b[0], 1./8./3.) == 2*howm*3);
            CHECK( num_cells( b[1], 2./8./2.) == 4*howm*2);
        }
        // Without X-point refine at the right x and the periodic y boundary
        for( dg::bc bcx : {dg::DIR, dg::PER})
        {
            INFO( "bcx "<<dg::bc2str( bcx));
            dg::EquidistXRefinement ref( 2, 1, 2, 2);
            auto b = check_refinedX( ref, 0., 0., bcx);
            CHECK( num_cells( b[0], 1./8./3.) == 2*3);
            CHECK( fabs( b[0][b[0].size()-7] - 0.75) < 1e-13);
            CHECK( num_cells( b[1], 2./8./2.) == 2*2*2);
        }
    }
    SECTION( "Exponential X Refinement")
    {
        // x is refined equidistantly, in y the cells next to the X-point are
        // halved add_y times, which leaves two cells of width h/2^add_y
        for( std::array<unsigned,3> p : {std::array<unsigned,3>{2,2,1},
                std::array<unsigned,3>{1,3,2}})
        {
            unsigned add_x = p[0], add_y = p[1], howm = p[2];
            INFO( "add_x "<<add_x<<" add_y "<<add_y<<" howmany "<<howm);
            dg::ExponentialXRefinement ref( add_x, add_y, howm, howm);
            auto b = check_refinedX( ref, 0.25, 0.25, dg::DIR);
            CHECK( num_cells( b[0], 1./8./(add_x+1)) == 2*howm*(add_x+1));
            CHECK( num_cells( b[1], 2./8./pow(2,add_y)) == 4*2);
        }
    }
}
