#pragma once
#include "dg/algorithm.h"
#include "tridiaginv.h" // lapack wrapper

/**
* @brief Functions for optimizing Contours
*/


namespace dg{
namespace mat{

// Problem: Newton method is attracted by saddle points. This we observe often i.e. Newton
// moves in the wrong diretion
// https://arxiv.org/pdf/1406.2572
// https://www.scientific.net/AMM.347-350.2586
//
template<class Func, class Jacobian, class InvHessian, class ContainerType>
unsigned newton( Func f, Jacobian jac, InvHessian invhess,
    ContainerType& x0, double tol = 1e-5, unsigned max_iter = 1000)
{
    ContainerType jj(x0), p(x0), test(x0);
    for ( unsigned i=0; i<max_iter; i++)
    {
        jac( x0, jj);
        invhess( x0, jj, p);
        double alpha = 1.;
        dg::blas1::axpby( -alpha, p, 1., x0);
        double err = sqrt(dg::blas1::dot( jj, jj));
        if (err < tol)
            return i;
    }
    return max_iter;
}

template<class Func, class Jacobian, class ContainerType>
unsigned levenberg_marquardt( Func fun, Jacobian jac,
    std::vector<double>& x0,
    ContainerType& copyable, // size of return of fun ( use somehow)
    double tol = 1e-8, unsigned max_iter = 1000)
{
    unsigned num_p = x0.size();
    auto x1 = x0, W = x0;
    auto rs(copyable), rs1(rs);
    std::vector<ContainerType> jacs(num_p, copyable);
    dg::SquareMatrix<double> HH(num_p, 0.), evHH( HH), evHH_T(HH), WW(HH);
    thrust::host_vector<double> evs( num_p), grad(num_p), gradbar(num_p),
        pk(num_p), pkbar(num_p);
    thrust::host_vector<double> work( 3*num_p-1);
    // init loop
    fun( x0, rs);
    jac( x0, jacs); // TODO FIX BUG HERE
    double delta = 0;
    for( unsigned p=0; p<num_p; p++)
        WW(p,p) = W[p] = dg::blas1::dot( jacs[p], jacs[p]);
    for( unsigned p=0; p<num_p; p++)
    {
        grad[p] = -dg::blas1::dot( jacs[p], rs);
        delta += grad[p]*grad[p]/W[p];
    }
    //std::cout << "Norm gT g/W^2 " << sqrt(delta)<<"\n";
    double f0 = dg::blas1::dot( rs, rs);
    delta = 0.25*f0/sqrt(delta);
    std::cout << "initial delta " << delta<<"\n";
    double normx0 = sqrt(dg::blas1::dot( x0, x0));
    for ( unsigned k=0; k<max_iter; k++)
    {
        // 1. Solve (J^T J + lambda W)p = -J^T r with lambda : ||p||_W leq Delta
        // W[p] = max{ W_{k-1}, H_pp}
        // In: jacs, rs; Out: pk, normpk, lambda
        for( unsigned l=0; l<num_p; l++)
        {
            for( unsigned j=l; j<num_p; j++)
                HH(j,l) = HH(l,j) = dg::blas1::dot( jacs[l], jacs[j]);
            WW(l,l) = std::max( WW(l,l), HH(l,l));
        }
        lapack::sygv( LAPACK_COL_MAJOR, 1, 'V', 'U', num_p, HH.data(), num_p, WW.data(), num_p, evs, work);
        evHH_T = HH;
        evHH = HH.transpose();
        //std::cout << "#########Iteration "<<k<<"\n";
        //std::cout << "Eigenvalues are \n";
        //for( unsigned p=0; p<num_p; p++)
        //    std::cout << "p EV "<<evs[p]<<"\n";
        //std::cout << "Weights are \n";
        //for( unsigned p=0; p<num_p; p++)
        //    std::cout << "W "<<WW(p,p)<<"\n";
        dg::blas2::gemv( evHH_T, grad, gradbar); // !! sygv gives A = WE_A Lambda E_A^TW
        double normp=0;
        auto target = [&]( double lambda)
        {
            // safeguard against 0 Eigenvalue!
            for( unsigned p=0; p<num_p; p++)
                pkbar[p]= gradbar[p]/(evs[p] +lambda == 0 ? 1e-16 : evs[p]+lambda);
            normp = sqrt(dg::blas1::dot( pkbar, pkbar));
            //std::cout << "Norm p "<<normp<<" delta "<<delta<<"\n";
            return 1./delta - 1./normp;
        };
        auto dtarget = [&]( double lambda)
        {
            // safeguard against 0 Eigenvalue!
            double dnorm =0;
            for( unsigned p=0; p<num_p; p++)
                dnorm += pkbar[p]*pkbar[p]/(evs[p]+lambda == 0 ? 1e-16 : evs[p]+lambda);
            return - dnorm/normp/normp/normp;
        };
        double lambda = 0;
        const double sigma = 1e-4; // tolerance for Newton algorithm
        const unsigned max_newton = 100;
        for( unsigned i=0; i<max_newton; i++) // Safeguard
        {
            double tt = target(lambda);
            //std::cout << "tt "<<tt<<" lambda "<<lambda<<"\n";
            // if p leq delta ( 1+sigma)
            if ( tt <= sigma/normp || i == max_newton-1) // tt can be negative
                break;
            lambda += -tt/dtarget(lambda);
        }
        dg::blas2::gemv(evHH, pkbar, pk);
        // target(lambda) updates grad and normp
        // 2. Check termination
        //std::cout << "Norm p "<<normp<<" normx0 "<<normx0<<"\n";
        //std::cout << "Real Norm p "<<sqrt(dg::blas1::dot( pk, pk))<<" normx0 "<<normx0<<"\n";
        if( sqrt(dg::blas1::dot( pk,pk)) <= tol*(normx0 + 1.))
        {
            dg::blas1::axpby( 1., pk , 1., x0);
            return k;
        }
        // 3. Compute Ratio rhok
        // don't overwrite (x0, rs) because we might reject
        dg::blas1::axpby( 1., pk , 1., x0, x1);
        //for( unsigned l=0; l<num_p; l++)
        //    std::cout << "x0 "<<l<<" "<<x0[l]<<"\n";
        //for( unsigned l=0; l<num_p; l++)
        //    std::cout << "x1 "<<l<<" "<<x1[l]<<"\n";
        //for( unsigned l=0; l<num_p; l++)
        //    std::cout << "pk "<<l<<" "<<pk[l]<<"\n";
        //for( unsigned l=0; l<num_p; l++)
        //    std::cout << "pkbar "<<l<<" "<<pkbar[l]<<"\n";
        //std::cout << "x0 "<<x0[0]<<"\n";
        //std::cout << "x1 "<<x1[0]<<"\n";
        fun( x1, rs1);
        double f1 = dg::blas1::dot( rs1, rs1);
        // compute (Jp)^2
        double jp = dg::blas2::dot( pkbar, evs, pkbar);
        //std::cout<< "f0 "<<f0<<" f1 "<<f1<<" "<<f0-f1<<"\n";
        double rhok = (1.-f1/f0)/( jp/f0 + 2*lambda*normp*normp/f0);
        //std::cout << "Ratio "<<rhok<<"\n";
        // Check
        ContainerType tmp(rs.size(),0.);
        for( unsigned i=0; i<num_p; i++)
            dg::blas1::axpby( pk[i], jacs[i], 1., tmp);
        //std::cout << "Real ratio "<<(f0-f1)/(-2*dg::blas1::dot( rs, tmp) - dg::blas1::dot( tmp, tmp))<<"\n";
        //std::cout << "Real diff  "<<(f0-f1)<<"\n";
        //std::cout << "Appr diff  "<<(-2*dg::blas1::dot( rs, tmp) - dg::blas1::dot( tmp, tmp))<<"\n";
        //std::cout << "Appr diff0 "<<(-2*dg::blas1::dot( rs, tmp))<<"\n";
        //std::cout << "Appr diff1 "<<(- dg::blas1::dot( tmp, tmp))<<"\n";
        if( rhok < 0.25)
            delta = 0.25*delta; // More has a more elaborate scheme
        else
        {
            // More has a bit different conditions
            if( rhok > 0.75 && lambda > 0)
                delta = 2*delta;
            // else delta remains unchanged
        }
        const double eta = 1e-4; // when step is rejected
        if ( rhok > eta)
        {
            x0 = x1;
            // update all quantities
            f0 = f1;
            using std::swap;
            swap( rs, rs1);
            jac( x0, jacs);
            normx0 = sqrt(dg::blas1::dot( x0, x0));
            for( unsigned l=0; l<num_p; l++)
            {
                grad[l] = -dg::blas1::dot( jacs[l], rs);
                //std::cout << "grad "<<l<<" "<<grad[l]<<"\n";
            }
        }
        // else step is rejected
    }
    return max_iter;
}



double f_alpha( double alphabar, double lambda = 1.)
{
    return 1-exp( -lambda*alphabar);
}
double finv_alpha( double alpha, double lambda = 1.)
{
    return - log( 1. - alpha)/lambda;
}
double df_alpha( double alphabar, double lambda = 1)
{
    return lambda*exp( -alphabar*lambda);
}
double ddf_alpha( double alphabar, double lambda = 1)
{
    return -lambda*lambda*exp( -alphabar*lambda);
}

std::pair<std::vector<thrust::complex<double>>,std::vector<thrust::complex<double>>>
    weights_and_nodes_talbot( unsigned N, const std::vector<double>& params)
{
    thrust::complex<double> I( 0,1);
    double h = M_PI/(double)N;
    unsigned n = N/2;
    std::vector<thrust::complex<double>> zk(N/2);
    std::vector<thrust::complex<double>> wk(N/2);
    double mu = params[0];
    double sigma = params[1];
    double nu = params[2];
    double alphabar = params[3];
    double alpha = f_alpha( alphabar);

    for( unsigned k=0; k<n; k++)
    {
        auto x = thrust::complex<double>( 2*h*k + h);
        zk[k] = thrust::complex<double>(N)*(-sigma + mu*x/tan(alpha*x) + nu*I*x);
        auto vk = thrust::complex<double>(N)*(mu/tan(alpha*x) - mu*x*alpha/sin(alpha*x)/sin(alpha*x) + nu*I);
        wk[k] = I*h/M_PI*vk;
    }
    return std::make_pair( zk, wk);
}

std::pair<std::vector<thrust::complex<double>>,std::vector<thrust::complex<double>>>
    jacobian_talbot( unsigned N, const std::vector<double>& params)
{
    thrust::complex<double> I( 0,1);
    double h = M_PI/(double)N;
    unsigned n = N/2;
    std::vector<thrust::complex<double>> dzk(4*n);
    std::vector<thrust::complex<double>> dwk(4*n);
    double mu = params[0];
    double alphabar = params[3];
    double alpha = f_alpha( alphabar);
    double dalpha = df_alpha( alphabar);
    for( unsigned k=0; k<n; k++)
    {
        auto Nc = thrust::complex<double>(N);
        auto x = thrust::complex<double>( 2*h*k + h);
        dzk[0*n+k] = Nc*( x/tan(alpha*x));
        dzk[1*n+k] = Nc*(-1. );
        dzk[2*n+k] = Nc*(I*x);
        dzk[3*n+k] = Nc*(-mu*x*x/sin(alpha*x)/sin(alpha*x));
        dzk[3*n+k] *= dalpha;
        dwk[0*n+k] = I*h/M_PI*Nc*(1./tan(alpha*x) - x*alpha/sin(alpha*x)/sin(alpha*x) );
        dwk[1*n+k] = 0;
        dwk[2*n+k] = I*h/M_PI*Nc*(I);
        dwk[3*n+k] = I*h/M_PI*Nc*(-2.*mu*x/sin(alpha*x)/sin(alpha*x) +
                             mu*x*x*alpha*2.*cos(alpha*x)/sin(alpha*x)/sin(alpha*x)/sin(alpha*x) );
        dwk[3*n+k] *= dalpha;

    }
    return std::make_pair( dzk, dwk);
}
std::pair<std::vector<thrust::complex<double>>,std::vector<thrust::complex<double>>>
    hessian_talbot( unsigned N, const std::vector<double>& params)
{
    thrust::complex<double> I( 0,1);
    double h = M_PI/(double)N;
    unsigned n = N/2;
    std::vector<thrust::complex<double>> ddzk(16*n, 0);
    std::vector<thrust::complex<double>> ddwk(16*n, 0);
    double mu = params[0];
    double alphabar = params[3];
    double alpha = f_alpha( alphabar);
    double dalpha = df_alpha( alphabar);
    double ddalpha = ddf_alpha( alphabar);
    if( alpha == 0)
    {
        std::cerr << "Warning alpha = 0 not allowed\n";
        alpha = 1e-6;
    }
    if( alpha > 1)
    {
        std::cerr << "Warning alpha > 1 not allowed\n";
        alpha = 1;
    }
    for( unsigned k=0; k<n; k++)
    {
        auto Nc = thrust::complex<double>(N);
        auto x = thrust::complex<double>( 2.*h*k + h);
        ddzk[(0*4+3)*n+k] = Nc*( -x*x/sin(alpha*x)/sin(alpha*x))*dalpha;
        ddzk[(3*4+0)*n+k] = Nc*( -x*x/sin(alpha*x)/sin(alpha*x))*dalpha;
        ddzk[(3*4+3)*n+k] = Nc*(mu*x*x*x*2.*cos(alpha*x)/sin(alpha*x)/sin(alpha*x)/sin(alpha*x));
        auto dzk = Nc*(-mu*x*x/sin(alpha*x)/sin(alpha*x));
        ddzk[(3*4+3)*n+k] = ddzk[(3*4+3)*n+k]*dalpha*dalpha + dzk*ddalpha;

        ddwk[(0*4+3)*n+k] = I*h/M_PI*Nc*(-2.*x/sin(alpha*x)/sin(alpha*x) +
                             x*x*alpha*2.*cos(alpha*x)/sin(alpha*x)/sin(alpha*x)/sin(alpha*x) )*dalpha;
        ddwk[(3*4+0)*n+k] = I*h/M_PI*Nc*(-2.*x/sin(alpha*x)/sin(alpha*x) +
                             x*x*alpha*2.*cos(alpha*x)/sin(alpha*x)/sin(alpha*x)/sin(alpha*x) )*dalpha;
        ddwk[(3*4+3)*n+k] = I*h/M_PI*Nc*(6.*mu*x*x*cos(alpha*x)/sin(alpha*x)/sin(alpha*x)/sin(alpha*x) -
                             mu*x*x*x*alpha*2.*(2.*cos(alpha*x)*cos(alpha*x)+1.)/sin(alpha*x)/sin(alpha*x)/sin(alpha*x)/sin(alpha*x));
        auto dwk = I*h/M_PI*Nc*(-2.*mu*x/sin(alpha*x)/sin(alpha*x) +
                             mu*x*x*alpha*2.*cos(alpha*x)/sin(alpha*x)/sin(alpha*x)/sin(alpha*x) );
        ddwk[(3*4+3)*n+k] = ddwk[(3*4+3)*n+k]*dalpha*dalpha + dwk*ddalpha;
    }
    return std::make_pair( ddzk, ddwk);
}

std::vector<double> weights_and_nodes2params( const
std::pair<std::vector<thrust::complex<double>>,std::vector<thrust::complex<double>>>& zkwk)
{
    const auto& zk = zkwk.first;
    const auto& wk = zkwk.second;
    std::vector<double> params( 4*zk.size());
    unsigned N = params.size()/2;
    for( unsigned i=0; i<zk.size(); i++)
    {
        params[2*i] = zk[i].real();
        params[2*i+1] = zk[i].imag();
    }
    for( unsigned i=0; i<zk.size(); i++)
    {
        params[N+2*i] = wk[i].real();
        params[N+2*i+1] = wk[i].imag();
    }
    return params;
}


std::pair<std::vector<thrust::complex<double>>,std::vector<thrust::complex<double>>>
    weights_and_nodes_identity( unsigned N, const std::vector<double>& params )
{
    if( N != params.size()/2)
        throw dg::Error(dg::Message(_ping_)<<"N "<<N<<" must match 0.5 params.size "<<params.size()/2<<"!");
    std::vector<thrust::complex<double>> zk(N/2);
    std::vector<thrust::complex<double>> wk(N/2);
    for( unsigned i=0; i<zk.size(); i++)
        zk[i] = thrust::complex<double>(params[2*i], params[2*i+1]);
    for( unsigned i=0; i<zk.size(); i++)
        wk[i] = thrust::complex<double>(params[N+2*i], params[N+2*i+1]);
    return std::make_pair( zk, wk);
}
std::pair<std::vector<thrust::complex<double>>,std::vector<thrust::complex<double>>>
    jacobian_identity( unsigned N, const std::vector<double>& params)
{
    unsigned n = N/2;
    std::vector<thrust::complex<double>> dzk(4*n*n, {0.});
    std::vector<thrust::complex<double>> dwk(4*n*n, {0.});
    for( unsigned k=0; k<n; k++)
    {
        dzk[(2*k+0)*n + k] = thrust::complex<double>(1,0);
        dzk[(2*k+1)*n + k] = thrust::complex<double>(0,1);
        dwk[(2*n+2*k+0)*n + k] = thrust::complex<double>(1,0);
        dwk[(2*n+2*k+1)*n + k] = thrust::complex<double>(0,1);
    }
    return std::make_pair( dzk, dwk);
}


/////////////////////////////////////////////////////////////////////////////

struct LeastSquaresCauchyError
{
    /**! @brief
     *
     * @param func must accept complex values as arguments
     */
    template<class Generator, class UnaryFunction>
    LeastSquaresCauchyError( unsigned N, Generator generate, UnaryFunction func, const
        std::vector<double>& rrs, const std::vector<double>& lls):
        m_N(N), m_func( func), m_generate(generate),
        m_rrs(rrs), m_lls(lls), m_exact( lls.size()*rrs.size()),
        m_func_rrs( N/2*rrs.size())
        {
            m_nl= lls.size();
            m_nr= rrs.size();
            for( unsigned i=0; i<m_nl; i++)
                for( unsigned j=0; j<m_nr; j++)
                    m_exact[i*m_nr+j] = (func( -thrust::complex<double>(m_lls[i]*m_rrs[j]))).real();
        }
    void result( const std::vector<double>& params, std::vector<double>& result)
    {
        auto pair = m_generate( m_N, params);
        const auto& zk = pair.first;
        const auto& wk = pair.second;
        unsigned n = zk.size();
        dg::blas1::copy( 0, result);
        for( unsigned k=0; k<n; k++)
            for( unsigned j=0; j<m_nr; j++)
                m_func_rrs[k*m_nr + j] = wk[k]*m_func( m_rrs[j]*zk[k]);
        for( int k=n-1; k>=0; k--)
            for( unsigned i=0; i<m_nl; i++)
                for( unsigned j=0; j<m_nr; j++)
                    result[i*m_nr+j]+= 2*(m_func_rrs[k*m_nr+j]/(-m_lls[i] - zk[k])).real();
    }
    void error( const std::vector<double>& params, std::vector<double>& err)
    {
        result( params, err);
        dg::blas1::axpby( -1., m_exact, 1., err);
    }

    void operator()( const std::vector<double>& params, std::vector<double>& res)
    {
        //dg::Timer t;
        //t.tic();
        result( params, res);
        //t.toc();
        //std::cout << "Computing result "<<t.diff()<<"\n";
        //t.tic();
        dg::blas1::axpby( -1., m_exact, 1., res);
        dg::blas1::pointwiseDot( res, res, res); // least squares converge better in r^4

        //t.toc();
        //std::cout << "Computing error  "<<t.diff()<<"\n";
    }
    private:
    unsigned m_N, m_nl, m_nr;
    std::function<thrust::complex<double>(thrust::complex<double>)> m_func;
    std::function<std::pair<std::vector<thrust::complex<double>>, std::vector<thrust::complex<double>>>(unsigned, const std::vector<double>&)> m_generate;
    std::vector<double> m_rrs, m_lls, m_exact;
    std::vector<thrust::complex<double>> m_func_rrs;
};

struct LeastSquaresCauchyJacobian
{
    /**! @brief
     *
     * @param func must accept complex values as arguments
     */
    template<class Generator, class GeneratorJac, class UnaryFunction, class UnaryFunctionD>
    LeastSquaresCauchyJacobian( unsigned N, Generator generate, GeneratorJac generateJac,
        UnaryFunction func, UnaryFunctionD dxlnfunc, const std::vector<double>& rrs, const std::vector<double>& lls):
        m_N(N), m_func( func), m_dxlnfunc(dxlnfunc),
        m_generate(generate), m_generateJac( generateJac),
        m_rrs(rrs), m_lls(lls), m_exact( lls.size()*rrs.size()),
        m_result( m_exact),
        m_func_rrs(N/2*rrs.size())
        {
            m_nl= lls.size();
            m_nr= rrs.size();
            for( unsigned i=0; i<m_nl; i++)
                for( unsigned j=0; j<m_nr; j++)
                    m_exact[i*m_nr+j] = (func( -thrust::complex<double>(m_lls[i]*m_rrs[j]))).real();
        }

    void operator()( const std::vector<double>& params, std::vector<std::vector<double>>& jac)
    {
        //dg::Timer t;
        //t.tic();
        auto pair = m_generate( m_N, params);
        const auto& zk = pair.first;
        const auto& wk = pair.second;
        auto Jacpair = m_generateJac( m_N, params);
        const auto& dzk = Jacpair.first;
        const auto& dwk = Jacpair.second;
        unsigned n = zk.size();
        //t.toc();
        //std::cout <<"Generating pairs took "<<t.diff()<<"\n";
        //t.tic();
        dg::blas1::copy( 0, m_result);
        for( unsigned k=0; k<n; k++)
            for( unsigned j=0; j<m_nr; j++)
                m_func_rrs[k*m_nr+j] = wk[k]*m_func( m_rrs[j]*zk[k]);
        for( int k=n-1; k>=0; k--)
            for( unsigned i=0; i<m_nl; i++)
                for( unsigned j=0; j<m_nr; j++)
                    m_result[i*m_nr+j]+= 2*(m_func_rrs[k*m_nr+j]/(-m_lls[i] - zk[k])).real();
        dg::blas1::axpby( -1., m_exact, 1., m_result);
        //t.toc();
        //std::cout <<"Computing result "<<t.diff()<<"\n";
        //t.tic();
        std::vector<thrust::complex<double>> tmp( m_func_rrs.size());
        dg::blas1::copy( 0, jac);
        for( unsigned p=0; p<params.size(); p++)
        {
            for( unsigned k=0; k<n; k++)
                for( unsigned j=0; j<m_nr; j++)
                    tmp[k*m_nr+j] =  dwk[p*n+k]/wk[k] +
                            m_rrs[j]*dzk[p*n+k]*m_dxlnfunc(m_rrs[j]*zk[k]);
            for( int k=n-1; k>=0; k--)
                for( unsigned i=0; i<m_nl; i++)
                    for( unsigned j=0; j<m_nr; j++)
                        jac[p][i*m_nr+j]+= 2*(m_func_rrs[k*m_nr+j]*(
                        dzk[p*n+k]/(-m_lls[i] - zk[k]) +
                        tmp[k*m_nr+j])/(-m_lls[i] - zk[k])).real();
            dg::blas1::pointwiseDot( 2., jac[p], m_result, 0., jac[p]);
        }
        //t.toc();
        //std::cout <<"Computing jacobian "<<t.diff()<<"\n";
    }
    private:
    unsigned m_N, m_nl, m_nr;
    std::function<thrust::complex<double>(thrust::complex<double>)> m_func, m_dxlnfunc;
    std::function<std::pair<std::vector<thrust::complex<double>>, std::vector<thrust::complex<double>>>(unsigned, const std::vector<double>&)> m_generate, m_generateJac;
    std::vector<double> m_rrs, m_lls, m_exact, m_result;
    std::vector<thrust::complex<double>> m_func_rrs;
};

struct CauchyHessian
{
    /**! @brief
     *
     * @param func must accept complex values as arguments
     */
    template<class Generator, class GeneratorJac, class GeneratorHessian,
        class UnaryFunction, class UnaryFunctionD, class UnaryFunctionDD>
    CauchyHessian( unsigned N, Generator generate, GeneratorJac generateJac,
        GeneratorHessian generateHess, UnaryFunction func, UnaryFunctionD dxlnfunc,
        UnaryFunctionDD dxxlnfunc, const std::vector<double>& rrs, const
        std::vector<double>& lls ):
        m_N(N), m_func( func), m_dxlnfunc(dxlnfunc), m_dxxlnfunc( dxxlnfunc),
        m_generate(generate), m_generateJac( generateJac), m_generateHess(generateHess),
        m_rrs(rrs), m_lls(lls), m_exact( lls.size()*rrs.size()),
        m_result( m_exact), m_hess(m_result),
        m_func_rrs(N/2*rrs.size())
        {
            m_nl= lls.size();
            m_nr= rrs.size();
            for( unsigned i=0; i<m_nl; i++)
                for( unsigned j=0; j<m_nr; j++)
                    m_exact[i*m_nr+j] = (func( -thrust::complex<double>(m_lls[i]*m_rrs[j]))).real();
        }

    void operator()( const std::vector<double>& params, const std::vector<double>& rhs, std::vector<double>& hessInvRhs)
    {
        auto pair = m_generate( m_N, params);
        const auto& zk = pair.first;
        const auto& wk = pair.second;
        auto Jacpair = m_generateJac( m_N, params);
        const auto& dzk = Jacpair.first;
        const auto& dwk = Jacpair.second;
        auto Hesspair = m_generateHess( m_N, params);
        const auto& ddzk = Hesspair.first;
        const auto& ddwk = Hesspair.second;
        unsigned n = zk.size();
        dg::blas1::copy( 0, m_result);
        assert( m_func_rrs.size() == n*m_nr);
        for( unsigned k=0; k<n; k++)
            for( unsigned j=0; j<m_nr; j++)
                m_func_rrs[k*m_nr+j] = wk[k]*m_func( m_rrs[j]*zk[k]);
        for( int k=n-1; k>=0; k--)
            for( unsigned i=0; i<m_nl; i++)
                for( unsigned j=0; j<m_nr; j++)
                    m_result[i*m_nr+j]+= 2*(m_func_rrs[k*m_nr+j]/(-m_lls[i] - zk[k])).real();
        dg::blas1::axpby( -1., m_exact, 1., m_result);
        double norm = dg::blas1::dot( m_result, m_result)/2.;
        std::vector<thrust::complex<double>> tmp( m_func_rrs.size());
        unsigned ps = params.size();
        m_jj.resize( ps);
        m_jac.resize( ps);
        for( unsigned p=0; p<ps; p++)
        {
            m_jac[p].resize( m_nl*m_nr);
            for( unsigned k=0; k<n; k++)
                for( unsigned j=0; j<m_nr; j++)
                    tmp[k*m_nr+j] = dwk[p*n+k]/wk[k] +
                            m_rrs[j]*dzk[p*n+k]*m_dxlnfunc(m_rrs[j]*zk[k]);
            dg::blas1::copy( 0, m_jac[p]);
            for( int k=n-1; k>=0; k--)
                for( unsigned i=0; i<m_nl; i++)
                    for( unsigned j=0; j<m_nr; j++)
                        m_jac[p][i*m_nr+j]+= 2*(m_func_rrs[k*m_nr+j]*(
                        dzk[p*n+k]/(-m_lls[i] - zk[k]) + tmp[k*m_nr+j]
                            )/(-m_lls[i] - zk[k])).real();
            m_jj[p] = dg::blas1::dot( m_result, m_jac[p]);
        }
        m_hh.resize( params.size()*params.size());
        std::vector<thrust::complex<double>> tmpw( n);
        std::vector<thrust::complex<double>> tmpp( m_func_rrs.size());
        std::vector<thrust::complex<double>> tmpq( m_func_rrs.size());
        for( unsigned p=0; p<ps; p++)
        for( unsigned q=p; q<ps; q++) // Hessian is symmetric
        {
            dg::blas1::copy( 0, m_hess);
            for( unsigned k=0; k<n; k++)
                for( unsigned j=0; j<m_nr; j++)
                {
                    tmpp[k*m_nr+j] = dwk[p*n+k]/wk[k] +
                            m_rrs[j]*dzk[p*n+k]*m_dxlnfunc(m_rrs[j]*zk[k]);
                    tmpq[k*m_nr+j] = dwk[q*n+k]/wk[k] +
                            m_rrs[j]*dzk[q*n+k]*m_dxlnfunc(m_rrs[j]*zk[k]);
                    tmp[k*m_nr+j] = m_rrs[j]*ddzk[(p*4+q)*n+k]*m_dxlnfunc(m_rrs[j]*zk[k])
                        +m_rrs[j]*m_rrs[j]*dzk[p*n+k]*dzk[q*n+k]*m_dxxlnfunc(m_rrs[j]*zk[k]);
                }
            for( unsigned k=0; k<n; k++)
                tmpw[k] = ddwk[(p*4+q)*n+k]/wk[k] -
                        dwk[p*n+k]*dwk[q*n+k]/wk[k]/wk[k];
            for( int k=n-1; k>=0; k--)
                for( unsigned i=0; i<m_nl; i++)
                    for( unsigned j=0; j<m_nr; j++)
                        m_hess[i*m_nr+j] += 2*( m_func_rrs[k*m_nr+j]*( (
                        tmpp[k*m_nr+j]+dzk[p*n+k]/( -m_lls[i] -zk[k]))*(
                        tmpq[k*m_nr+j]+dzk[q*n+k]/( -m_lls[i] -zk[k])) +
                        tmpw[k] + tmp[k*m_nr+j]
                        +ddzk[(p*4+q)*n+k]/(-m_lls[i]-zk[k]) +
                        dzk[p*n+k]*dzk[q*n+k]/(-m_lls[i]-zk[k])/(-m_lls[i]-zk[k]))/(-m_lls[i]-
                        zk[k] )).real();
            m_hh[p*ps+q] = dg::blas1::dot( m_jac[p],m_jac[q]);
            m_hh[p*ps+q] += dg::blas1::dot( m_result, m_hess);
            m_hh[p*ps+q] = m_hh[p*ps+q]/norm - m_jj[p]*m_jj[q]/norm/norm;
        }
        for( unsigned p=0; p<ps; p++)
            for( unsigned q=0; q<p; q++) // Hessian is symmetric
                m_hh[p*ps+q] = m_hh[q*ps+p];
        std::vector<unsigned> pivots;
        dg::Operator<double> HH (m_hh);
        dg::create::lu_pivot( HH,  pivots);
        hessInvRhs = rhs;
        dg::create::lu_solve( HH, pivots, hessInvRhs);
    }
    std::vector<double> last_hessian( ) const{ return m_hh;}
    private:
    unsigned m_N, m_nl, m_nr;
    std::function<thrust::complex<double>(thrust::complex<double>)> m_func, m_dxlnfunc, m_dxxlnfunc;
    std::function<std::pair<std::vector<thrust::complex<double>>,
        std::vector<thrust::complex<double>>>(unsigned, const std::vector<double>&)>
            m_generate, m_generateJac, m_generateHess;
    std::vector<double> m_rrs, m_lls, m_exact, m_result, m_hess, m_jj, m_hh;
    std::vector<std::vector<double>> m_jac;
    std::vector<thrust::complex<double>> m_func_rrs;
};


std::vector<double> generate_range( double min, double max, unsigned per_order = 20)
{
    unsigned orders = unsigned (log10(max) - log10(min));
    if ( orders == 0)
        orders = 1;
    std::vector<double> range(orders*per_order);
    unsigned N = range.size();
    double h = (log10(max) - log10(min))/double(N-1);
    for( unsigned i=0; i<N; i++)
        range[i] = pow( 10.0, log10(min) + i*h);
    return range;
}
///////////////////////////////////////////////////////////////////////////////////////

// Not such a great idea:
// Mainly because of the lack of the accuracy of our blas1::dot functions
template<class ContainerType, class UnaryFunc>
void result_talbot( unsigned N, double rrs, double lls,
    const std::vector<ContainerType>& ps,
    UnaryFunc func,
    ContainerType& result)
{
    dg::blas1::subroutine( [N,rrs,lls,func]DG_DEVICE(
        double mu, double sigma, double nu, double alphabar, double& result)
    {
        thrust::complex<double> I( 0,1);
        double h = M_PI/(double)N;
        unsigned n = N/2;
        result = 0;
        double alpha = 1.-exp(-alphabar);
        for( int k=n-1; k>=0; k--)
        {
            double x = 2*h*k + h;
            double tanx = tan(alpha*x);
            double sinx = sin(alpha*x);
            thrust::complex<double> zk( N*(-sigma + mu*x/tanx), nu*x*N);
            thrust::complex<double> wk(-h/M_PI*N*nu, h*N/M_PI*(mu/tanx - mu*x*alpha/sinx/sinx));
            result+= 2*(wk*func(rrs*zk)/(-lls - zk)).real();
        }
    }, ps[0], ps[1], ps[2], ps[3], result);
}

template<class ContainerType, class UnaryFunc>
void error_talbot( unsigned N, const std::vector<double>& rrs, const std::vector<double>& lls,
    const std::vector<ContainerType>& ps,
    UnaryFunc f,
    ContainerType& error)
{
    dg::blas1::copy( 0, error);
    ContainerType tmp( error);
    for( unsigned i=0; i<lls.size();i++)
    for( unsigned j=0; j<rrs.size();j++)
    {
        result_talbot( N, rrs[j], lls[i], ps, f, tmp);
        dg::blas1::plus( tmp, -f(-thrust::complex<double>(lls[i]*rrs[j])).real());
        dg::blas1::pointwiseDot( 1., tmp, tmp, 1., error);
    }
}

} //namespace mat
} //namespace dg

