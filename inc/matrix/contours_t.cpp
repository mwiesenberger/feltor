#include <iostream>
#include "functors.h"

#include "contours.h"


int main()
{
    unsigned n = 7;
    std::vector<double> params = {0.5017,0.6122,0.2645,dg::mat::finv_alpha(0.6407)};
    auto pair = dg::mat::weights_and_nodes_talbot( 2*n, params);
    const auto& zk = pair.first;
    const auto& wk = pair.second;
    for( unsigned i=0; i<zk.size(); i++)
        std::cout << zk[i]<<" "<<wk[i]<<"\n";
    std::cout << "################################\n\n";
    std::vector<double> rrs = {1};
    std::vector<double> lls = {0,1,2,10};
    auto func = dg::mat::GyrolagK<thrust::complex<double>>(0,1);
    auto dxlnfunc = dg::mat::DLnGyrolagK<thrust::complex<double>>(0,1);
    auto dxxlnfunc = dg::mat::DDLnGyrolagK<thrust::complex<double>>(0,1);
    dg::mat::LeastSquaresCauchyError cauchy( 2*n, dg::mat::weights_and_nodes_talbot,
        func, rrs, lls);
    std::vector<double> results( lls.size()*rrs.size());
    cauchy.result( params, results);
    for( unsigned i=0; i<results.size(); i++)
    {
        std::cout << "i "<<i<<" "<<results[i]<<"\n";
        double result;
        result_talbot( 2*n, rrs[0], lls[i], params, func, result);
        std::cout << "i "<<i<<" "<<result<<"\n";

    }
    cauchy( params, results);
    std::cout << "Cauchy error "<<dg::blas1::dot( results, results)<<"\n";
    double error;
    error_talbot( 2*n, rrs, lls, params, func, error);
    std::cout << "GPU func error "<<error<<"\n";
    std::cout << "################################\n\n";
    rrs = dg::mat::generate_range( 0.0625, 12.5);
    lls = dg::mat::generate_range( 21, 23401947);
    dg::mat::LeastSquaresCauchyError cauchyTCV( 2*n, dg::mat::weights_and_nodes_talbot,
        func, rrs, lls);
    dg::mat::LeastSquaresCauchyJacobian jacTCV( 2*n, dg::mat::weights_and_nodes_talbot,
        dg::mat::jacobian_talbot,
        func, dxlnfunc, rrs, lls);
    dg::mat::CauchyHessian hessTCV( 2*n, dg::mat::weights_and_nodes_talbot,
        dg::mat::jacobian_talbot, dg::mat::hessian_talbot,
        func, dxlnfunc, dxxlnfunc, rrs, lls);
    std::cout << "problem size is "<<lls.size()*rrs.size()<<"\n";
    results.resize( lls.size()*rrs.size());
    auto results_eps = results;
    cauchyTCV( params, results);
    std::cout << "Cauchy error "<<dg::blas1::dot( results, results)<<"\n";
    std::vector<std::vector<double>> jac(4, {results});
    jacTCV( params, jac);
    std::cout << "Cauchy Jacobian\n";
    for( unsigned i=0; i<params.size(); i++)
    {
        double gradC = 2*dg::blas1::dot( jac[i], results);
        std::cout << gradC<<std::endl;
    }
    std::cout << "Numerical Jacobian\n";
    for( unsigned i=0; i<params.size(); i++)
    {
        std::vector<double> eps( params.size(), 0);
        double tol = 1e-6;
        eps[i] = tol;
        dg::blas1::axpby( 1., params, 1., eps);
        cauchyTCV( eps, results_eps);
        dg::blas1::axpby( 1, results_eps, -1, results, results_eps);
        double gradN = 2*dg::blas1::dot( results_eps, results)/tol;
        std::cout << gradN<<std::endl;
    }
    std::cout << "\n";
    auto zkwk = dg::mat::weights_and_nodes_talbot( 2*n, params);
    auto paramsI = dg::mat::weights_and_nodes2params( zkwk);
    dg::mat::LeastSquaresCauchyError IcauchyTCV( 2*n, dg::mat::weights_and_nodes_identity,
        func, rrs, lls);
    dg::mat::LeastSquaresCauchyJacobian IjacTCV( 2*n, dg::mat::weights_and_nodes_identity,
        dg::mat::jacobian_identity,
        func, dxlnfunc, rrs, lls);
    std::cout << "Cauchy Talbot   error "<<dg::blas1::dot( results, results)<<"\n";
    IcauchyTCV( paramsI, results);
    std::cout << "Cauchy Identity error "<<dg::blas1::dot( results, results)<<"\n";
    std::vector<std::vector<double>> jacI(paramsI.size(), {results});
    IjacTCV( paramsI, jacI);
    std::cout << "Cauchy Identity Jacobian\n";
    for( unsigned i=0; i<paramsI.size(); i++)
    {
        double gradC = 2*dg::blas1::dot( jacI[i], results);
        std::cout << gradC<<std::endl;
    }
    std::cout << "Numerical Identity Jacobian\n";
    for( unsigned i=0; i<paramsI.size(); i++)
    {
        std::vector<double> eps( paramsI.size(), 0);
        double tol = 1e-10;
        eps[i] = tol;
        dg::blas1::axpby( 1., paramsI, 1., eps);
        IcauchyTCV( eps, results_eps);
        dg::blas1::axpby( 1, results_eps, -1, results, results_eps);
        double gradN = 2*dg::blas1::dot( results_eps, results)/tol;
        std::cout << gradN<<std::endl;
    }
    std::cout << "\n";
    dg::Timer t;
    t.tic();
    for( unsigned i=0; i<1000; i++)
        cauchyTCV( params, results);
    t.toc();
    std::cout << "One cauchy eval took (s) "<<t.diff()/(double)1000<<"\n";
    t.tic();
    for( unsigned i=0; i<100; i++)
        jacTCV( params, jac);
    t.toc();
    std::cout << "One jacobian eval took (s) "<<t.diff()/(double)100<<"\n";
    std::vector<double> tmp(params);
    t.tic();
    for( unsigned i=0; i<10; i++)
        hessTCV( params, tmp, tmp);
    t.toc();
    std::cout << "One Hessian eval took (s) "<<t.diff()/(double)10<<"\n";

    unsigned steps = levenberg_marquardt( cauchyTCV, jacTCV, params, results, 1e-4, 1000);
    std::cout << "Found optimum in "<<steps<<" steps\n";
    cauchyTCV( params, results);
    std::cout << "Cauchy error "<<dg::blas1::dot( results, results)<<"\n";
    return 0;
}
