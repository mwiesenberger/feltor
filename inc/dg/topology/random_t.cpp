#include <iostream>

#include "random.h"

#include "catch2/catch_all.hpp"

TEST_CASE( "Uniform real distribution")
{
    dg::Grid2d g2d( 0, 1, 0, 1, 3, 10, 20);
    double a = -1, b = 2;
    auto vec = dg::create::uniform_real_distribution( a, b, g2d);
    REQUIRE( vec.size() == g2d.size());
    SECTION( "Range")
    {
        for( unsigned u=0; u<vec.size(); u++)
        {
            INFO( "Index "<<u);
            CHECK( a <= vec[u]);
            CHECK( vec[u] < b);
        }
    }
    SECTION( "Sequence of the default minstd_rand engine")
    {
        thrust::minstd_rand rng;
        thrust::uniform_real_distribution<double> dist( a, b);
        for( unsigned u=0; u<vec.size(); u++)
        {
            INFO( "Index "<<u);
            CHECK( vec[u] == dist( rng));
        }
    }
    SECTION( "Single elements with discard")
    {
        // The MPI version relies on this
        thrust::uniform_real_distribution<double> dist( a, b);
        for( unsigned u : {0u, 1u, 7u, 100u, (unsigned)vec.size()-1})
        {
            thrust::minstd_rand rng;
            rng.discard( u);
            INFO( "Index "<<u);
            CHECK( vec[u] == dist( rng));
        }
    }
    SECTION( "Value type and dimensions")
    {
        dg::RealGrid3d<float> g3d( 0, 1, 0, 1, 0, 1, 2, 3, 4, 5);
        auto vecf = dg::create::uniform_real_distribution( 0, 1, g3d);
        static_assert( std::is_same_v<decltype(vecf), thrust::host_vector<float>>);
        CHECK( vecf.size() == g3d.size());
        dg::Grid1d g1d( 0, 1, 2, 7);
        CHECK( dg::create::uniform_real_distribution( 0, 1, g1d).size() == g1d.size());
    }
}
