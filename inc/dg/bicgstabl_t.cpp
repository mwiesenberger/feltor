#include <iostream>
#include <iomanip>

#ifdef WITH_MPI
#include <mpi.h>
#include "backend/mpi_init.h"
#endif

#include "pcg.h"
#include "bicgstabl.h"
#include "andersonacc.h"
#include "lgmres.h"
#include "elliptic.h"
#include "chebyshev.h"
#include "catch2/catch_all.hpp"

const double lx = 2.*M_PI;
const double ly = 2.*M_PI;


double fct(double x, double y){ return sin(y)*sin(x);}
double laplace_fct( double x, double y) { return 2*sin(y)*sin(x);}
double initial( double, double) {return sin(0);}

TEST_CASE( "Solvers")
{
#ifdef WITH_MPI
    int rank,size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm comm = dg::mpi_cart_create( MPI_COMM_WORLD, {0,0}, {1,1});
#endif
    const unsigned n=4, Nx=36, Ny = 48;
    //global relative error in L2 norm is O(h^P)
    //more N means less iterations for same error
    INFO( "Computing on the Grid " <<n<<" x "<<Nx<<" x "<<Ny);
    dg::x::Grid2d grid( 0, lx, 0, ly,n, Nx, Ny, dg::DIR, dg::PER
#ifdef WITH_MPI
        , comm
#endif
    );
    unsigned max_iter = n*n*Nx*Ny;
    const dg::x::HVec solution = dg::evaluate ( fct, grid);
    // create volume and inverse volume on previously defined grid
    const dg::x::HVec w2d = dg::create::weights( grid);
    // Evaluate right hand side and solution on the grid
    const dg::x::HVec b = dg::evaluate ( laplace_fct, grid);
    double lmin = 2.0, lmax = 2*grid.size(); //Eigenvalues of Laplace

    dg::x::HVec x = dg::evaluate( initial, grid);
    const dg::x::HVec& copyable_vector = x;
    dg::x::HVec error( solution);
    // Create unnormalized Laplacian
    dg::Elliptic<dg::x::CartesianGrid2d, dg::x::HMatrix, dg::x::HVec> A( grid, dg::forward);

    double res;
    res = sqrt(dg::blas2::dot(w2d , solution));
    INFO( "L2 Norm of Solution is        " << res);
    SECTION("PCG")
    {
        INFO("PCG SOLVER");
        //! [pcg]
        // allocate memory in conjugate gradient
        dg::PCG pcg( copyable_vector, max_iter);

        // Solve
        unsigned num_iter = pcg.solve( A, x, b, 1., w2d, 1e-6);

        INFO( "Number of pcg iterations "<< num_iter);
        dg::blas1::axpby( 1.,x,-1.,error);
        res = sqrt(dg::blas2::dot(w2d , error));
        INFO( "L2 Norm of Error is           " << res);
        CHECK( res < 1e-6);
        //! [pcg]
    }
    SECTION( "cheby" )
    {
        // Chebyshev converges slowly here because the condition number
        // lmax/lmin is large; we only check that the error decreases
        INFO("CHEBYSHEV SOLVER");
        const double res0 = res;
        dg::ChebyshevIteration cheby( x);
        unsigned num_iter =200;
        cheby.solve( A, x, b, lmin, 2*lmax, num_iter);
        INFO("After "<<num_iter<<" Chebyshev iterations");
        dg::blas1::axpby( 1.,x,-1.,error);
        res = sqrt(dg::blas2::dot(w2d , error));
        INFO( "L2 Norm of Error is           " << res);
        CHECK( res < 0.5*res0);
    }
    SECTION( "P cheby" )
    {
        // Chebyshev converges slowly here because the condition number
        // lmax/lmin is large; we only check that the error decreases
        INFO("PRECONDITIONED CHEBYSHEV SOLVER");
        const double res0 = res;
        dg::ChebyshevIteration cheby( x);
        unsigned num_iter =200;
        cheby.solve( A, x, b, A.precond(), lmin, lmax, num_iter);
        INFO( "After "<<num_iter<<" Chebyshev iterations");
        dg::blas1::axpby( 1.,x,-1.,error);
        res = sqrt(dg::blas2::dot(w2d , error));
        INFO( "L2 Norm of Error is           " << res);
        CHECK( res < 0.5*res0);
    }
    SECTION( "bicgstabl" )
    {
        INFO("BICGSTABl SOLVER");
        //! [bicgstabl]
        dg::BICGSTABl bicg( x, 100, 2);
        unsigned num_iter = bicg.solve( A, x, b, A.precond(), A.weights(), 1e-6);
        INFO( "After "<<num_iter<<" BICGSTABl iterations we have");
        dg::blas1::axpby( 1.,x,-1.,error);
        res = sqrt(dg::blas2::dot(w2d , error));
        INFO( "L2 Norm of Error is           " << res);
        CHECK( res < 1e-6);
        //! [bicgstabl]
    }
    SECTION( "lgmres" )
    {
        INFO("LGMRES SOLVER");
        //! [lgmres]
        dg::LGMRES lgmres( x, 30, 4, 10000);
        unsigned num_iter = lgmres.solve( A, x, b, A.precond(), A.weights(), 1e-6);
        INFO( "After "<<num_iter<<" LGMRES iterations we have");
        dg::blas1::axpby( 1.,x,-1.,error);
        res = sqrt(dg::blas2::dot(w2d , error));
        INFO( "L2 Norm of Error is           " << res);
        CHECK( res < 1e-6);
        //! [lgmres]
    }
    SECTION( "Andersonacc" )
    {
        INFO("ANDERSONACC");
        //! [andersonacc]
        dg::AndersonAcceleration acc( copyable_vector, 10);
        const double eps = 1e-6;
        double damping = 1e-3;
        unsigned restart = 10;
        unsigned num_iter = acc.solve( A, x, b, w2d, eps, eps, max_iter,
                damping, restart, false);
        INFO( "Number of iterations "<< num_iter);
        dg::blas1::axpby( 1.,x,-1.,error);
        res = sqrt(dg::blas2::dot(w2d , error));
        INFO( "L2 Norm of Error is           " << res);
        CHECK( res < 1e-6);
        //! [andersonacc]
    }


}

TEST_CASE( "BICGSTABl")
{
    std::vector<double> x( 4, 0.);
    CHECK_THROWS_AS( dg::BICGSTABl<std::vector<double>>( x, 10, 0), dg::Error);
}

TEST_CASE( "Chebyshev")
{
    // A diagonal matrix with Eigenvalues in [1,3] and the identity as weights
    std::vector<double> A = {1., 1.5, 2., 3.}, b = {1., 1., 1., 1.}, x( 4, 0.);
    dg::ChebyshevIteration<std::vector<double>> cheby( x);
    const double lmin = 1., lmax = 3., theta = 2.;
    SECTION( "Convergence")
    {
        cheby.solve( A, x, b, lmin, lmax, 40);
        for( unsigned u=0; u<4; u++)
            CHECK( fabs( x[u] - 1./A[u]) < 1e-12);
        // The preconditioner P = 2 scales the spectrum of PA to [2,6]
        x.assign( 4, 0.);
        cheby.solve( A, x, b, 2., 2*lmin, 2*lmax, 40);
        for( unsigned u=0; u<4; u++)
            CHECK( fabs( x[u] - 1./A[u]) < 1e-12);
    }
    SECTION( "One iteration")
    {
        // One iteration from x=0 is x_1 = P b / theta with or without x_is_zero
        std::vector<double> x0( 4, 0.), x1( 4, NAN), y0( 4, 0.), y1( 4, NAN);
        cheby.solve( A, x0, b, lmin, lmax, 1, false);
        cheby.solve( A, x1, b, lmin, lmax, 1, true);
        cheby.solve( A, y0, b, 1., lmin, lmax, 1, false);
        cheby.solve( A, y1, b, 1., lmin, lmax, 1, true);
        for( unsigned u=0; u<4; u++)
        {
            CHECK( x0[u] == 1./theta);
            CHECK( x1[u] == 1./theta);
            CHECK( y0[u] == 1./theta);
            CHECK( y1[u] == 1./theta);
        }
    }
    SECTION( "ModifiedChebyshevPreconditioner")
    {
        // The polynomial approximates the inverse of A with an error of
        // O(q^degree) with q = (sqrt(lmax/lmin)-1)/(sqrt(lmax/lmin)+1) = 0.27
        dg::ModifiedChebyshevPreconditioner<std::vector<double>, std::vector<double>>
            pre( A, x, lmin, lmax, 20);
        dg::blas2::symv( pre, b, x);
        for( unsigned u=0; u<4; u++)
            CHECK( fabs( x[u] - 1./A[u]) < 1e-10);
    }
    SECTION( "x_is_zero")
    {
        // x_is_zero ignores the content of x (also NaN)
        std::vector<double> x0( 4, 0.), x1( 4, NAN), y0( 4, 0.), y1( 4, NAN);
        cheby.solve( A, x0, b, lmin, lmax, 5, false);
        cheby.solve( A, x1, b, lmin, lmax, 5, true);
        cheby.solve( A, y0, b, 1., lmin, lmax, 5, false);
        cheby.solve( A, y1, b, 1., lmin, lmax, 5, true);
        for( unsigned u=0; u<4; u++)
        {
            CHECK( fabs( x1[u] - x0[u]) < 1e-15);
            CHECK( fabs( y0[u] - x0[u]) < 1e-15);
            CHECK( fabs( y1[u] - x0[u]) < 1e-15);
        }
    }
}
