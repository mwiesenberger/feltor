#include <iostream>
#include <thrust/complex.h>
#include <complex>
#include "functors.h"

#include "catch2/catch_all.hpp"

TEST_CASE( "phi functions")
{
    // phi_k(x) = sum_j x^j/(j+k)! is accurate everywhere, in particular close
    // to 0 and around |x| = 2, where the evaluation switches from the series
    // to the recursion (reference values in quadruple precision)
    std::vector<double> limits = {1,0.5,1./6.,1./24.};
    std::vector<double> xs = {-10., -2.01, -1.99, -0.0101, 1e-10, 1e-5, 1e-3,
        0.0099, 0.5, 1., 1.99, 2.01, 10.};
    std::vector<std::array<double,4>> refs = {
        {9.99954600070237515e-02, 9.00004539992976248e-02, 4.09999546000702375e-02, 1.25666712066596429e-02}, // -10
        {4.30851405637410496e-01, 2.83158504658004758e-01, 1.07881340966166798e-01, 2.92464306967661067e-02}, // -2.01
        {4.33821394248480475e-01, 2.84511862186693230e-01, 1.08285496388596367e-01, 2.93372714965177385e-02}, // -1.99
        {9.94966958824029571e-01, 4.98320908511923673e-01, 1.66246681987755142e-01, 4.15826414763885509e-02}, // -0.0101
        {1.00000000005000000e+00, 5.00000000016666667e-01, 1.66666666670833333e-01, 4.16666666675000000e-02}, // 1e-10
        {1.00000500001666671e+00, 5.00001666670833342e-01, 1.66667083334166668e-01, 4.16667500001388891e-02}, // 1e-05
        {1.00050016670834167e+00, 5.00166708341668056e-01, 1.66708341668055754e-01, 4.16750013890873264e-02}, // 0.001
        {1.00496637550930694e+00, 5.01654091849185504e-01, 1.67079984766212472e-01, 4.17493029844248149e-02}, // 0.0099
        {1.29744254140025629e+00, 5.94885082800512587e-01, 1.89770165601025175e-01, 4.62069978687170162e-02}, // 0.5
        {1.71828182845904524e+00, 7.18281828459045235e-01, 2.18281828459045235e-01, 5.16151617923785687e-02}, // 1
        {3.17363505643696808e+00, 1.09227892283264728e+00, 2.97627599413390594e-01, 6.58095139430773509e-02}, // 1.99
        {3.21558076981054397e+00, 1.10227898995549463e+00, 2.99641288535071986e-01, 6.61565282927389720e-02}, // 2.01
        {2.20254657948067165e+03, 2.20154657948067165e+02, 2.19654657948067165e+01, 2.17987991281400499e+00}}; // 10

    std::vector<std::function<double(double)>> funcs = {dg::mat::phi1<double>, dg::mat::phi2<double>, dg::mat::phi3<double>, dg::mat::phi4<double>};
    unsigned u = GENERATE( 0, 1, 2, 3);
    INFO( "Phi"<<u+1);
    INFO( "Ana "<<limits[u] << "\nNum "<<funcs[u](0));
    CHECK( limits[u] == funcs[u](0));
    for( unsigned i=0; i<xs.size(); i++)
    {
        INFO( "x "<<xs[i]<<" Num "<<funcs[u](xs[i])<<" Ana "<<refs[i][u]);
        CHECK( fabs( funcs[u](xs[i]) - refs[i][u]) < 2e-15*refs[i][u]);
    }
}

TEST_CASE( "Gyrolag functions")
{
    double a = 0.5;

    dg::mat::GyrolagK gyro0( 0, a), gyro1( 1, a), gyro2( 2, a);
    CHECK ( (gyro0( 1.) - exp(-0.5)) < 1e-12);
    CHECK ( (gyro1( 1.) - 2.*exp(-0.5)) < 1e-12);
    CHECK ( (gyro2( 1.) - 2.*exp(-0.5)) < 1e-12);

    dg::mat::DGyrolagK dgyro0( 0, a), dgyro1( 1, a), dgyro2( 2, a);
    CHECK ( (dgyro0( 1.) - (-0.303265)) < 1e-6);
    CHECK ( (dgyro1( 1.) - 0.151633) < 1e-6);
    CHECK ( (dgyro2( 1.) - 0.113724) < 1e-6);


}

TEST_CASE( "Complex Gyrolag functions")
{
    // On the real axis the complex functions equal the real ones
    double a = 0.5, x = 1.3;
    for( unsigned n : {0u, 1u, 2u, 3u, 4u, 7u})
    {
        INFO( "n "<<n);
        double ref = dg::mat::GyrolagK<double>( n, a)( x);
        double dref = dg::mat::DGyrolagK<double>( n, a)( x);
        thrust::complex<double> t = dg::mat::GyrolagK<thrust::complex<double>>( n, a)( x);
        std::complex<double> s = dg::mat::GyrolagK<std::complex<double>>( n, a)( x);
        CHECK( fabs( t.real() - ref) < 1e-15);
        CHECK( t.imag() == 0);
        CHECK( fabs( s.real() - ref) < 1e-15);
        CHECK( s.imag() == 0);
        thrust::complex<double> dt = dg::mat::DGyrolagK<thrust::complex<double>>( n, a)( x);
        std::complex<double> ds = dg::mat::DGyrolagK<std::complex<double>>( n, a)( x);
        CHECK( fabs( dt.real() - dref) < 1e-15);
        CHECK( fabs( ds.real() - dref) < 1e-15);
    }
}
