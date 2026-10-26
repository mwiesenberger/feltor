#include <iostream>
#include "average.h"
#include "stencil.h"
#include "filter.h"
#include "../blas2.h"

#include "catch2/catch_all.hpp"

// Index and sign of point i of a window stencil in a grid with N points
static unsigned window_index( int i, unsigned N, dg::bc bc, double& sign)
{
    if( i >= 0 && i < (int)N)
        return i;
    if( bc == dg::PER)
        return i < 0 ? i + N : i - N;
    if( bc == dg::DIR)
        sign = -1;
    return i < 0 ? -i-1 : 2*N-1-i;
}
// Does the slope limiter flatten cell c of a linear increasing function?
static bool limiter_flattens( unsigned c, unsigned N, dg::bc bc)
{
    if( c != 0 && c != N-1)
        return false;
    // PER: the neighbour across the boundary is on the other end
    // NEU: the mirrored neighbour has the same mean
    // DIR: the mirrored neighbour has the negative mean, which is
    // consistent with the slope at the left (x0 = 0) but not the right
    if( bc == dg::DIR)
        return c == N-1;
    return true;
}

TEST_CASE( "Stencil")
{
    std::vector<dg::bc> bcs = {dg::DIR, dg::NEU, dg::PER};
    SECTION( "Window stencil")
    {
        dg::Grid2d g2d( 0,1, 0,2, 3, 5, 3);
        const unsigned Nx = g2d.gx().size(), Ny = g2d.gy().size();
        dg::HVec x( g2d.size()), y(x), z(x);
        // integer values make all sums exact
        for( unsigned i=0; i<x.size(); i++)
            x[i] = i;
        for( auto bc : bcs)
        {
            INFO( "Test "<<dg::bc2str( bc)<<" boundary:");
            auto stencil = dg::create::window_stencil( {3,3}, g2d, bc, bc);
            dg::blas2::symv( stencil, x, y);
            dg::blas2::stencil( dg::CSRSymvFilter(), (dg::IHMatrix)stencil, x, z);
            for( unsigned i=0; i<Ny; i++)
            for( unsigned k=0; k<Nx; k++)
            {
                // Sum over the 3x3 window; boundary points are duplicated,
                // with a minus sign for DIR, or wrapped for PER
                double sol = 0;
                for( int a=-1; a<=1; a++)
                for( int b=-1; b<=1; b++)
                {
                    double sy = 1, sx = 1;
                    unsigned ii = window_index( (int)i+a, Ny, bc, sy);
                    unsigned kk = window_index( (int)k+b, Nx, bc, sx);
                    sol += sy*sx*x[ii*Nx+kk];
                }
                INFO( "Row "<<i<<" column "<<k);
                CHECK( y[i*Nx+k] == sol);
                CHECK( z[i*Nx+k] == sol);
            }
        }
    }

    SECTION( "dg slope limiter stencil")
    {
        dg::Grid1d g1d( 0,1, 3, 3);
        const dg::HVec x = dg::evaluate( dg::cooX1d, g1d);
        for( auto bc : bcs)
        {
            INFO( "Test "<<dg::bc2str( bc)<<" boundary");
            auto stencil = dg::create::limiter_stencil( g1d, bc);
            // The zero coefficient of each cell has 3n entries, the others 0
            for( unsigned i=0; i<g1d.size(); i++)
                CHECK( stencil.row_offsets()[i+1] - stencil.row_offsets()[i]
                        == ( i%g1d.n() == 0 ? 3*(int)g1d.n() : 0));
            dg::HVec y(x);
            dg::blas2::stencil( dg::CSRSlopeLimiter<double>(), stencil, x, y);
            for( unsigned c=0; c<g1d.N(); c++)
            {
                // A linear function is kept unless the slope from the
                // neighbour cell has the opposite sign; then the cell is
                // flattened to its mean
                bool flat = limiter_flattens( c, g1d.N(), bc);
                double mean = g1d.x0() + (c+0.5)*g1d.h();
                for( unsigned j=0; j<g1d.n(); j++)
                {
                    unsigned idx = c*g1d.n()+j;
                    INFO( "Cell "<<c<<" point "<<j);
                    CHECK( fabs( y[idx] - (flat ? mean : x[idx])) < 1e-14);
                }
            }
        }
    }
    SECTION( "Test DIR boundary")
    {
        dg::Grid2d g2d( 0,1, 0,2, 3, 5, 3);
        auto stencil = dg::create::limiter_stencil( dg::coo3d::y, g2d, dg::DIR);
        auto x = dg::evaluate( [](double x, double y){return sin(x)*y;}, g2d);
        auto y = x;
        dg::blas2::stencil( dg::CSRSlopeLimiter<double>(), stencil, x, y);
        // The limiter acts along y only: the last cell is flattened
        // to its mean, the others are kept
        const unsigned Nx = g2d.gx().size();
        for( unsigned c=0; c<g2d.Ny(); c++)
        for( unsigned j=0; j<g2d.ny(); j++)
        for( unsigned k=0; k<Nx; k++)
        {
            unsigned idx = (c*g2d.ny()+j)*Nx+k;
            double mean = sin( g2d.gx().abscissas()[k])*(g2d.y0()+(c+0.5)*g2d.hy());
            bool flat = limiter_flattens( c, g2d.Ny(), dg::DIR);
            INFO( "Cell "<<c<<" point "<<j<<" column "<<k);
            CHECK( fabs( y[idx] - (flat ? mean : x[idx])) < 1e-14);
        }
    }


    SECTION( "Transpose")
    {
        dg::Grid2d g2d( 0,1, 0,2, 3, 5, 3);
        auto x = dg::evaluate( dg::one, g2d), y(x);
        auto xx(x), yy(x), zz(x);
        dg::Grid2d g2dT( 0,2, 0,1, 3, 3, 5);
        auto stencilX = dg::create::limiter_stencil( dg::coo3d::x, g2dT, dg::DIR);
        dg::blas2::stencil( dg::CSRSlopeLimiter<double>(), stencilX, xx, zz);
        //for( unsigned i=0; i<g2d.gy().size(); i++)
        //    std::cout<< zz[i]<<" ";
        //std::cout << "\n";
        // transpose
        unsigned nx = g2d.shape(0), ny = g2d.shape(1);
        dg::blas2::parallel_for( [ny,nx]DG_DEVICE( unsigned k, const
                    double* ii, double* oo)
        {
            unsigned i = k/nx, j =  k%nx;
            oo[j*ny+i] = ii[i*nx+j];
        }, nx*ny, zz, yy);
        dg::blas1::axpby( 1., yy, -1., y);
        double err = dg::blas1::dot( y, 1.);
        INFO( "error between transposed application and app is: " <<err<<" (0)");
        CHECK( err == 0);
    }
}
