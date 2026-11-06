#include <iostream>
#include "dg/algorithm.h"
#include "exp_runge_kutta.h"

#include "catch2/catch_all.hpp"

// integrate Dot y = -100 y +sin(t)
const double matrix = -1;

static void rhs( double t, double, double& yp)
{
    //yp = 0.;
    yp = sin(t);
}

struct MatrixFunction
{
    template<class UnaryOp>
    void operator() ( UnaryOp f, double x, double& y) const{
        y = f(matrix)*x;
    }
};

static double solution( double t, double y0){
    return y0*exp(matrix*t) + (exp(matrix*t) - matrix*sin(t) -cos(t))/(1.+matrix*matrix);
    //return y0*exp(matrix*t);
}

// Dot y = matrix y + y^2: the nonlinearity depends on y, so the stage
// coefficients a_ij of the tableau matter
static void rhs_nonlinear( double, double y, double& yp)
{
    yp = y*y;
}
static double solution_nonlinear( double t, double y0){
    return matrix/( (matrix/y0 + 1.)*exp(-matrix*t) - 1.);
}


TEST_CASE( "Exp Runge Kutta")
{
    double t_start = 0., y_start = 1.;
    double t_end = M_PI/2.;
    const unsigned N = 40;
    const double dt = (t_end - t_start)/(double)N;

    double u = y_start;
    double u1(u), sol = y_start*exp(matrix*t_end);
    dg::mat::ExponentialStep<double> rk( u);
    dg::SinglestepTimeloop<double>( rk, MatrixFunction(),
            dt).integrate( t_start, u, t_end, u1);
    dg::blas1::axpby( 1., sol , -1., u1);
    double error = sqrt( dg::blas1::dot( u1, u1));
    INFO( "Norm of error in "<<std::setw(24) <<"ExponentialStep"<<"\t"<<error);
    CHECK( error < 1e-15);

    auto name = GENERATE( as<std::string>{},
        "Euler",
        "Midpoint-2-2",
        "Runge-Kutta-4-4",
        "Hochbruck-3-3-4"
    );
    auto b = dg::mat::create::func_tableau<double>(name);
    std::vector<unsigned> NTs = {20,40};
    std::vector<double> err(NTs.size());
    SECTION( "Nonlinearity independent of y")
    {
        for( unsigned k = 0; k<NTs.size(); k++)
        {
            double u = solution(t_start, y_start);
            double u1(u), sol = solution(t_end, y_start);
            dg::mat::ExponentialERKStep<double> rk( name, u);
            auto mat = MatrixFunction();
            const double dt = (t_end - t_start)/(double)NTs[k];
            dg::SinglestepTimeloop<double>( rk, std::tie(rhs,mat),
                dt).integrate( t_start, u, t_end, u1);
            dg::blas1::axpby( 1., sol , -1., u1);
            err[k] = sqrt( dg::blas1::dot( u1,u1));
        }
    }
    SECTION( "Nonlinearity depends on y")
    {
        for( unsigned k = 0; k<NTs.size(); k++)
        {
            double u = 0.5, u1(u), sol = solution_nonlinear( t_end, u);
            dg::mat::ExponentialERKStep<double> rk( name, u);
            auto mat = MatrixFunction();
            const double dt = (t_end - t_start)/(double)NTs[k];
            dg::SinglestepTimeloop<double>( rk, std::tie(rhs_nonlinear,mat),
                dt).integrate( t_start, u, t_end, u1);
            err[k] = fabs( sol - u1);
        }
    }
    double order = log( err[0]/err[1])/log( (double)NTs[1]/(double)NTs[0]);
    INFO("Norm of error in "<<std::setw(24) <<name<<"\t"<<err[1]<<" order "
        <<order<<" expected "<<b.order());
    // Hochbruck-3-3-4 is an exponential Rosenbrock method: its order 4 needs
    // the matrix to be the Jacobian of the whole right hand side at the
    // current step (and an autonomous problem), which is not the case here
    if( name == "Hochbruck-3-3-4")
        CHECK( (2  - order)/double(2) < 0.014);
    else
        CHECK( (b.order()  - order)/double(b.order()) < 0.014);
}
