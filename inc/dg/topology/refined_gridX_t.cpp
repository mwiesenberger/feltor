#include <iostream>

#include "../blas.h"
#include "transform.h"
#include "refined_gridX.h"

#include "catch2/catch_all.hpp"

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
}
