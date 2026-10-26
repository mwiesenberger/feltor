#include <iostream>
#include <iomanip>
#include <limits>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "blas.h"
#include "functors.h"

#include "catch2/catch_test_macros.hpp"
#include "catch2/generators/catch_generators.hpp"
#include "catch2/generators/catch_generators_range.hpp"

struct Expression{
   DG_DEVICE
   void operator() ( double& v, double w, double param){
       v = param*v*v + w;
   }
};


template<class Vector>
auto result( const Vector& vec) { return vec[0];}
#ifdef WITH_MPI
template<class Vector>
auto result( const dg::MPI_Vector<Vector>& vec) { return vec.data()[0];}
#endif

TEST_CASE( "blas")
{
#ifdef WITH_MPI
    int rank,size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm comm = MPI_COMM_WORLD;
#else
    int size = 1;
#endif
    //"This program tests the many different possibilities to call blas1 and blas2 functions."

    double x1 = 1., x2 = 2., x3 = 3.;
    std::vector<double> vec1(3, x3), vec2;
    vec1[0] = 10., vec1[1] = 20., vec1[2] = 30;
    vec2 = vec1;
    std::array<float, 3>  arr1{2,3,4}, arr2(arr1);
#ifdef WITH_MPI
    dg::MHVec hvec1 ( dg::HVec(vec1.begin(), vec1.end()), comm);
    dg::MDVec dvec1 ( dg::DVec(vec1.begin(), vec1.end()), comm);
#else
    dg::HVec  hvec1 ( vec1.begin(), vec1.end());
    dg::DVec  dvec1 ( vec1.begin(), vec1.end());
#endif
    std::vector<dg::x::DVec  > arrdvec1( 3, dvec1);

    SECTION( "trivial parallel functions")
    {
        dg::blas1::axpby( 2., x1, 3., x2);
        INFO( "Scalar addition");
        CHECK(x2 == 8.);

        INFO( "Reduction");
        CHECK( dg::blas1::reduce( arr1, 0, thrust::maximum<float>()) == 4);

        dg::blas1::axpby( 2., arr1, 3., vec2);
        INFO( "Recursive Vec Scalar addition");
        CHECK( vec2[0] == 34.);

        dg::blas1::axpby( 2., vec1, 3., arr2);
        INFO( "Recursive Arr Scalar addition");
        CHECK( arr2[0] == 26.);

        INFO( "Recursive DVec reduction");
        CHECK(dg::blas1::reduce( dvec1, 0, thrust::maximum<double>()) == 30);

        dg::blas1::copy( 2., arrdvec1);
        INFO( "Recursive DVec Copy Scalar to");
        CHECK(result(arrdvec1[0]) == 2);
        CHECK(result(arrdvec1[1]) == 2);

        //dg::blas1::axpby( 2., vec1 , 3, arrdvec1);
        for( unsigned i=0; i<3; i++)
            dg::blas1::axpby( 2., vec1[i], 3, arrdvec1[i]);
        INFO( "Recursive Scalar/Vetor addition");
        CHECK( result(arrdvec1[0]) == 26);
        CHECK( result(arrdvec1[1]) == 46.);

        // test the examples in the documentation
        // dg::blas1::subroutine( []__host__ __device__(double& v){ v+=1.;}, dvec1);
        dg::blas1::plus( dvec1, 1);
        std::array<dg::x::DVec, 3> array_v{ dvec1, dvec1, dvec1},
            array_w(array_v);
        INFO( result(dvec1)<< " "<<result(array_w[2]));
        dg::blas1::subroutine( Expression(), dvec1, array_w[2], 3);
        INFO( "Example in documentation          ");
        CHECK( result(dvec1) == 374);
        //dg::blas1::subroutine( Expression(), array_v, array_w, array_p);
        //INFO( "Example in documentation          "<< (array_v[0][0] == 132 && array_v[1][1] == 903)<<std::endl;
    }
    SECTION( "Test DOT functions:")
    {
        //Repeat result on arrdvec1
        dg::blas1::copy( 2., arrdvec1);
        for( unsigned i=0; i<3; i++)
            dg::blas1::axpby( 2., vec1[i], 3, arrdvec1[i]);

        std::array<double, 3> array_p{ 1,2,3};
        double result = dg::blas1::dot( 1., array_p);
        INFO( "blas1 dot recursive Scalar          ");
        CHECK(result == 6);

        result = dg::blas1::dot( 1., arrdvec1 )/(double)size;
        INFO( "blas1 dot recursive Scalar Vector   ");
        CHECK(result == 414);

        result = dg::blas2::dot( 1., 4., arrdvec1 )/(double)size;
        INFO( "blas2 dot recursive Scalar Vector   ");
        CHECK(result == 1656);

        result = dg::blas2::dot( 1., arrdvec1, arrdvec1 )/(double)size;
        double result1 = dg::blas1::dot( arrdvec1, arrdvec1)/(double)size;
        INFO( "blas1/2 dot recursive Scalar Vector ");
        CHECK(result == result1);
    }
    SECTION( "SYMV functions")
    {
        //Repeat result on arrdvec1
        dg::blas1::copy( 2., arrdvec1);
        for( unsigned i=0; i<3; i++)
            dg::blas1::axpby( 2., vec1[i], 3, arrdvec1[i]);

        // Repeat result on dvec1
        dg::blas1::plus( dvec1, 1);
        std::array<dg::x::DVec, 3> array_v{ dvec1, dvec1, dvec1},
            array_w(array_v);
        dg::blas1::subroutine( Expression(), dvec1, array_w[2], 3);

        dg::blas2::symv( 2., arrdvec1, arrdvec1);
        INFO( "symv Scalar times Vector          ");
        CHECK( result(arrdvec1[0]) == 52);
        std::array<std::vector<dg::x::DVec>,1> recursive{ arrdvec1};
        dg::blas2::symv( arrdvec1[0], recursive, recursive);
        INFO( "symv deep Recursion               ");
        CHECK( result(recursive[0][0]) == 52*52);
        dg::blas2::symv( 0.5, dg::asDenseMatrix(dg::asPointers(arrdvec1)),
                std::array<double,3>({0.1, 10, 1000}), 0.001, dvec1);
        INFO( "symv as DenseMatrix               ");
        CHECK( result(dvec1) == 66462.974);
        struct Functor{
            void operator()( const dg::x::DVec& in, dg::x::DVec&  out){
                dg::blas1::axpby( 1., in, 0., out);
            }
        };
        Functor f;
        dg::blas2::symv( f, arrdvec1[0], dvec1);
        INFO( "symv with functor ");
        CHECK( result(dvec1) == 52);
    }
    //Check compiler error:
    //struct NoFunctor{
    //
    //};
    //NoFunctor nof;
    //dg::blas2::symv( nof, arrdvec1[0], dvec1);
    SECTION( "std::map")
    {
        std::map< std::string, dg::x::DVec> testmap{ { "a", dvec1}, {"b", dvec1}};
        std::map< std::string, dg::x::DVec> testmap2( testmap);
        std::map< std::string, dg::x::DVec> testmap3{ { "c", dvec1}, {"b", dvec1}};
        CHECK_THROWS_AS(
            dg::blas1::axpby( 2., testmap, 3., testmap3),
            dg::Error);
        dg::blas1::axpby( 1., testmap, 1., testmap2);
        INFO( "axpby "<<result(testmap2["b"]) );
        CHECK( result(testmap2["a"]) == 20);
        CHECK( result(testmap2["b"]) == 20);
    }

}

// Distribute the same local vector to all ranks
template<class Vector>
static auto make_x_vector( const thrust::host_vector<dg::get_value_type<Vector>>& local)
{
#ifdef WITH_MPI
    return dg::MPI_Vector<Vector>( dg::construct<Vector>( local), MPI_COMM_WORLD);
#else
    return dg::construct<Vector>( local);
#endif
}

TEST_CASE( "Exact dot product")
{
    // The sum is exact before it is rounded: 2^100 + 1 - 2^100 = 1
    // Sizes test the remainder of the vectorized loop (multiples of 8) and
    // the OpenMP version (size >= 100)
    int ranks = 1;
#ifdef WITH_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &ranks);
#endif
    unsigned size = GENERATE( range( 1u, 21u), values( {99u, 100u, 101u, 1003u}));
    INFO( "Local size "<<size<<" on "<<ranks<<" ranks");
    const double big = ldexp( 1., 100);
    thrust::host_vector<double> h_v( size), h_w( size, 1.);
    unsigned plus = 0, minus = 0, ones = 0;
    for( unsigned i=0; i<size; i++)
    {
        if( i%3 == 0)
            h_v[i] = big, plus++;
        else if( i%3 == 1)
            h_v[i] = 1., ones++;
        else
            h_v[i] = -big, minus++;
    }
    // A single rounding of the exact result
    const double sol = double(plus-minus)*ranks*big + double(ones*ranks);
    auto v = make_x_vector<dg::DVec>( h_v);
    auto w = make_x_vector<dg::DVec>( h_w);
    SECTION( "Two vectors")
    {
        CHECK( dg::blas1::dot( v, w) == sol);
        CHECK( dg::blas1::dot( w, v) == sol);
    }
    SECTION( "Vector and scalar")
    {
        CHECK( dg::blas1::dot( v, 1.) == sol);
    }
    SECTION( "Three vectors")
    {
        CHECK( dg::blas2::dot( v, w, w) == sol);
        CHECK( dg::blas2::dot( w, v, w) == sol);
    }
    SECTION( "Negative sum")
    {
        dg::blas1::scal( v, -1.);
        CHECK( dg::blas1::dot( v, w) == -sol);
        CHECK( dg::blas2::dot( v, w, w) == -sol);
    }
    SECTION( "Float")
    {
        thrust::host_vector<float> h_fv( h_v.begin(), h_v.end());
        auto fv = make_x_vector<dg::fDVec>( h_fv);
        CHECK( dg::blas1::dot( fv, 1.f) == float(sol));
    }
}

TEST_CASE( "Dot of subnormal numbers")
{
    // Subnormal sums are exact in double precision as well
    int ranks = 1;
#ifdef WITH_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &ranks);
#endif
    // Sizes test the serial and the OpenMP version (size >= 100)
    unsigned size = GENERATE( 5u, 1003u);
    INFO( "Local size "<<size<<" on "<<ranks<<" ranks");
    const double dmin = std::numeric_limits<double>::min();
    const double tiny = std::numeric_limits<double>::denorm_min();
    thrust::host_vector<double> h_v( size, 0.), h_w( size, 1.);
    h_v[0] = tiny, h_v[1] = dmin/4., h_v[2] = dmin, h_v[3] = -dmin;
    const double sol = ranks*(tiny + dmin/4.);
    REQUIRE( sol < dmin); // the result is subnormal
    auto v = make_x_vector<dg::DVec>( h_v);
    auto w = make_x_vector<dg::DVec>( h_w);
    double dot2 = dg::blas1::dot( v, w), dot3 = dg::blas2::dot( v, w, w);
    INFO( std::hexfloat<<"Solution "<<sol<<" blas1 "<<dot2<<" blas2 "<<dot3);
    CHECK( dot2 == sol);
    CHECK( dot3 == sol);
}

TEST_CASE( "Dot throws on NaN or Inf")
{
    // Sizes test the serial and the OpenMP version (size >= 100)
    // With MPI only even ranks hold the non-finite value but all ranks throw
    int rank = 0;
#ifdef WITH_MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
    unsigned size = GENERATE( 10u, 1003u);
    double nonfinite = GENERATE( std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::infinity());
    INFO( "Size "<<size<<" with "<<nonfinite<<" on rank "<<rank);
    thrust::host_vector<double> h_v( size, 1.);
    if( rank%2 == 0)
        h_v[size/2] = nonfinite;
    auto v = make_x_vector<dg::DVec>( h_v);
    bool hasnan = dg::blas1::reduce( v, false,
            thrust::logical_or<bool>(), dg::ISNFINITE<double>());
    CHECK( hasnan);
    CHECK_THROWS_AS( dg::blas1::dot( v, v), dg::Error);
    CHECK_THROWS_AS( dg::blas2::dot( v, 1., v), dg::Error);
}
#ifndef WITH_MPI
TEST_CASE( "Blas2 documentation")
{
    SECTION( "dot")
    {
        //! [dot]
        dg::DVec two( 100,2), three(100,3);
        int result = dg::blas2::dot(two, 42., three);
        CHECK( result == 25200); //100*(2*42*3)
        //! [dot]
    }
    SECTION( "parallel_for")
    {
        //! [parallel_for]
        // Compute forward difference of vector
        std::vector<double> vec1 = {11,12,13}, vec2(3);
        dg::blas2::parallel_for(
            []DG_DEVICE( unsigned i, const double* x, double* y, int mod)
            {
                y[i] = x[(i+1)%mod] - x[i];
            }, 3, vec1, vec2, 3);
        CHECK( vec2 == std::vector{ 1.,1.,-2.});
        //! [parallel_for]
        //! [parallel_transpose]
        // Compute transpose of vector
        unsigned nx = 3, ny = 3;
        std::vector<double> vec = {11,12,13, 21,22,23, 31,32,33}, vecT(9);
        dg::blas2::parallel_for( [nx,ny]DG_DEVICE( unsigned k, const
                    double* ii, double* oo)
            {
                unsigned i = k/nx, j =  k%nx;
                oo[j*ny+i] = ii[i*nx+j];
            }, nx*ny, vec, vecT);
        CHECK( vecT == std::vector<double>{ 11,21,31, 12,22,32, 13,23,33});
        //! [parallel_transpose]
    }
}
#endif
