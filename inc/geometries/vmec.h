#pragma once
#include <string>
#include <vector>
#include "dg/algorithm.h"

namespace dg
{
namespace geo
{
namespace vmec
{
///@cond
namespace detail
{
struct VmecFunctorRTP
{
    VmecFunctorRTP() {}
    VmecFunctorRTP( const std::vector<double>& cos, // can be empty, access cos[ u * K + k], i.e. modes lie contiguous
                 const std::vector<double>& sin, // can be empty
                 const std::vector<int>& exponents, // for r
                 const std::vector<int>& mk, // poloidal modes for cos and sin
                 const std::vector<int>& nk) : // toroidal modes for cos and sin
                 m_cos(cos), m_sin(sin), m_exponents(exponents),
                 m_mk(mk), m_nk(nk), m_prev( {1e300,1e300,1e300,1e300}){
                     assert( m_mk.size() == m_nk.size());
     }

    double operator() (double r, double t, double p)
    {
        if( m_prev[2] == p && m_prev[1] == t && m_prev[0] == r)
            return m_prev[3];
        double result;
        unsigned K = m_mk.size();
        if( not m_cos.empty())
        {
            for( unsigned u=0; u<m_exponents.size(); u++)
                for( unsigned k=0; k<K; k++)
                {
                    result += m_cos[u*K + k]*pow( r, m_exponents[u])*cos( m_mk[k]*t - m_nk[k]*p);
                }
        }
        if( not m_sin.empty())
        {
            for( unsigned u=0; u<m_exponents.size(); u++)
                for( unsigned k=0; k<K; k++)
                {
                    result += m_sin[u*K + k]*pow( r, m_exponents[u])*sin( m_mk[k]*t - m_nk[k]*p);
                }
        }
        m_prev[0] = r, m_prev[1] = t, m_prev[2] = p, m_prev[3] = result;
        return result;
    }
    VmecFunctorRTP dr() const
    {
        std::vector<double> new_exp, new_cos, new_sin;
        unsigned K = m_mk.size();
        for( unsigned u=0; u<m_exponents.size(); u++)
            if( m_exponents[u] != 0)
            {
                new_exp.append( m_exponents[u] -1);
                for( unsigned k=0; k<K; k++)
                {
                    if( not m_cos.empty())
                        new_cos.append( m_exponent[u]*m_cos[u*K+k]);
                    if( not m_sin.empty())
                        new_sin.append( m_exponent[u]*m_sin[u*K+k]);
                }
            }
        return VmecFunctorRTP( new_cos, new_sin, new_exp, m_mk, m_nk);
    }
    VmecFunctorRTP dt() const
    {
        std::vector<double> new_cos, new_sin, new_mk, new_nk;
        unsigned K = m_mk.size();
        for( unsigned k=0; k<K; k++)
        {
            if( m_mk[k] != 0)
            {
                new_mk.append(m_mk[k]);
                new_nk.append(m_nk[k]);
            }
        }
        for( unsigned u=0; u<m_exponents.size(); u++)
            for( unsigned k=0; k<K; k++)
            {
                if( not m_cos.empty())
                    if( m_mk[k] != 0)
                    {
                        new_sin.append( -m_mk[k]*m_cos[u*K+k]);
                    }
                if( not m_sin.empty())
                    if( m_mk[k] != 0)
                        new_cos.append( m_mk[k]*m_sin[u*K+k]);
            }
        return VmecFunctorRTP( new_cos, new_sin, m_exponents, new_mk, new_nk);
    }
    VmecFunctorRTP dp() const
    {
        std::vector<double> new_cos, new_sin, new_mk, new_nk;
        unsigned K = m_mk.size();
        for( unsigned k=0; k<K; k++)
        {
            if( m_nk[k] != 0)
            {
                new_mk.append(m_mk[k]);
                new_nk.append(m_nk[k]);
            }
        }
        for( unsigned u=0; u<m_exponents.size(); u++)
            for( unsigned k=0; k<K; k++)
            {
                if( not m_cos.empty())
                    if( m_mk[k] != 0)
                    {
                        new_sin.append( m_nk[k]*m_cos[u*K+k]);
                    }
                if( not m_sin.empty())
                    if( m_mk[k] != 0)
                        new_cos.append( -m_nk[k]*m_sin[u*K+k]);
            }
        return VmecFunctorRTP( new_cos, new_sin, m_exponents, new_mk, new_nk);
    }

    private:
    std::vector<double> m_cos; // can be empty
    std::vector<double> m_sin; // can be empty
    std::vector<int> m_exponents; // for r
    std::vector<int> m_mk, m_nk; // for cos and sin
    mutable std::array<double,4> m_prev;
};

struct ConvertRZP2RTP
{
    ConvertRZP2RTP( const std::array<VmecFunctorRTP,2>& R0) :
        m_R0(R0[0]), m_Z0(R0[1]),
        m_R0p(R0[0].dp()), m_Z0p(R0[1].dp()),
        m_R0pp(R0[0].dp().dp()), m_Z0pp(R0[1].dp().dp()){}

    std::array<double,3> operator()( double R, double Z, double P) const
    {
        if( m_prev[2] == p && m_prev[1] == t && m_prev[0] == r)
            return m_prev_out;

        double R0 = m_R0(0,0,P);
        double Z0 = m_Z0(0,0,P);
        double Rbar = R - R0, Zbar = Z - Z0;
        double rad = sqrt( Rbar*Rbar + Zbar*Zbar);
        double the;
        if ( rad == 0)
            rad = 1e-20; // prevent nan
        if( Zbar > 0)
            the = arccos( Rbar / rad);
        else
            the = 2.*M_PI - arccos( Rbar / rad);
        m_prev[0] = R, m_prev[1] = Z, m_prev[2] = P;
        m_prev_out = {rad,the,P};
        return m_prev_out;
    }

    double rp( double r, double t, double p) const
    {
        return cos(t)*m_R0p(0,0,P) + sin(t)*m_Z0p(0,0,P);
    }
    double tp( double r, double t, double p) const
    {
        return (sin(t)*m_R0p(0,0,P) - cos(t)*m_Z0p(0,0,P))/r;
    }
    double rpp( double r, double t, double p) const
    {
        return cos(t)*m_R0pp(0,0,P) + sin(t)*m_Z0pp(0,0,P);
    }
    double tpp( double r, double t, double p) const
    {
        return (sin(t)*m_R0pp(0,0,P) - cos(t)*m_Z0pp(0,0,P))/r;
    }

    private:

    detail::VmecFunctorRTP m_R0, m_Z0, m_R0p, m_Z0p, m_R0pp, m_Z0pp;

    mutable std::array<double,3> m_prev, m_prev_out;


};

struct VmecFunctorsRTP
{
    VmecFunctorRTP v;
    VmecFunctorRTP vr, vt, vp;
    VmecFunctorRTP vrr, vrt, vrp, vtt, vtp, vpp;
};
VmecFunctorsRTP make_FunctorsRTP( const VmecFunctor& v)
{
    return {v,
     v.dr(), v.dt(), v.dp(),
     v.dr().dr(), v.dr().dt(), v.dr().dp(), v.dt().dt(), v.dt().dp(), v.dp()};
}


struct VmecFunctor : public aCylindricalFunctor<VmecFunctor>
{
    VmecFunctor( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto p = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return m_v->v(r,t,p);
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
struct VmecFunctorR : public aCylindricalFunctor<VmecFunctorR>
{
    VmecFunctorR( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return cos(t) * m_v->vr(r,t,p) - sin(t) * m_v->vt(r,t,p) / r;
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
struct VmecFunctorZ : public aCylindricalFunctor<VmecFunctorZ>
{
    VmecFunctorZ( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return sin(t) * m_v->vr(r,t,p) + cos(t) * m_v->vt(r,t,p) / r;
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
struct VmecFunctorP : public aCylindricalFunctor<VmecFunctorP>
{
    VmecFunctorP( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return m_v->vr(r,t,p) * m_c->rp(r,t,p) + m_v->vt(r,t,p) * m_c->tp(r,t,p) + m_v->vp(r,t,p);
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
struct VmecFunctorRR : public aCylindricalFunctor<VmecFunctorRR>
{
    VmecFunctorRR( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return (sin(t)*sin(t)*m_v->vtt(r,t,p) + sin(2*t)*m_v->vt(r,t,p))/r/r
          + (sin(t)*sin(t)*m_v->vr(r,t,p) - 2*sin(t)*cos(t)*m_v->vrt(r,t,p))/r + cos(t)*cos(t)*m_v->vrr(r,t,p);
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
struct VmecFunctorRZ : public aCylindricalFunctor<VmecFunctorRZ>
{
    VmecFunctorRZ( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return (-0.5*sin(2*t)*m_v->vtt(r,t,p) - cos(2*t)*m_v->vt(r,t,p))/r/r
          + (-0.5*sin(2*t)*m_v->vr(r,t,p) + cos(2*t)*m_v->vrt(r,t,p))/r + 0.5*sin(2*t)*m_v->vrr(r,t,p);
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
struct VmecFunctorRP : public aCylindricalFunctor<VmecFunctorRP>
{
    VmecFunctorRP( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return ( m_c->rp(r,t,p) * ( sin(t)*m_v->vt(r,t,p) + r*(r*cos(t)*m_v->vrr(r,t,p) - sin(t)*m_v->vrt(r,t,p)))
                -m_c->tp(r,t,p) * ( sin(t)*m_v->vtt(r,t,p) + r*sin(t)*m_v->vr(r,t,p) + cos(t)*m_v->vt(r,t,p) - r*cos(t)*m_v->vrt(r,t,p))
           + r*(r*cos(t)*m_v->vrp(r,t,p) - sin(t)*m_v->vtp(r,t,p)) ) / r /r;
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
struct VmecFunctorZZ : public aCylindricalFunctor<VmecFunctorZZ>
{
    VmecFunctorZZ( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return ( r*(cos(t)*cos(t) * m_v->vtt(r,t,p) + r*cos(t)*cos(t)*m_v->vr(r,t,p) + sin(t)*(r*sin(t)*m_v->vrr(r,t,p) + (r+1)*cos(t)*m_v->vrt(r,t,p))) - (r+1)*sin(t)*cos(t)*m_v->vt(r,t,p)  ) / r / r ;
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
struct VmecFunctorZP : public aCylindricalFunctor<VmecFunctorZP>
{
    VmecFunctorZP( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return ( m_c->rp(r,t,p) * ( -cos(t)*m_v->vt(r,t,p) + r*(r*sin(t)*m_v->vrr(r,t,p) + cos(t)*m_v->vrt(r,t,p)))
                +m_c->tp(r,t,p) * ( cos(t)*m_v->vtt(r,t,p) + r*cos(t)*m_v->vr(r,t,p) - sin(t)*m_v->vt(r,t,p) + r*sin(t)*m_v->vrt(r,t,p))
           + r*(r*sin(t)*m_v->vrp(r,t,p) + cos(t)*m_v->vtp(r,t,p)) ) / r /r;
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};

struct VmecFunctorPP : public aCylindricalFunctor<VmecFunctorPP>
{
    VmecFunctorPP( const std::shared_ptr<detail::VmecFunctorsRTP>& v, const std::shared_ptr<detail::ConvertRZP2RTP>& convert) :
     m_v(v), m_c(convert)
     {}

    double do_compute(double R, double Z, double P)
    {
        auto pp = *m_c(R,Z,P);
        double r = pp[0], t = p[1], p = p[2];
        return ( r*r*m_c->rpp(r,t,p)*m_v->vr(r,t,p)
                   + m_c->rp(r,t,p)*(2*r*r*m_v->vrp(r,t,p) - m_c->tp(r,t,p)*(m_v->vt(r,t,p) - 2*r*m_v->vrt(r,t,p)))
            + r*r*m_c->rp(r,t,p)*m_c->rp(r,t,p)*m_v->vrr(r,t,p) + r*r*m_v->vpp(r,t,p) - r*m_c->tp(r,t,p)*m_c->tp(r,t,p)*m_v->vr(r,t,p) + r*m_c->tpp(r,t,p)*m_v->vt(r,t,p) + m_c->tp(r,t,p)*m_c->rp(r,t,p)*m_v->vt(r,t,p) + 2*r*m_c->tp(r,t,p)*m_v->vtp(r,t,p) + m_c->tp(r,t,p)*m_c->tp(r,t,p)*m_v->vtt(r,t,p)  ) / r /r;
    }
    private:
    std::shared_ptr<detail::VmecFunctorsRTP> m_v;
    std::shared_ptr<detail::ConvertRZP2RTP> m_c;
};
}
///@endcond

struct ClebschParameters
{
    unsigned nfp; //!< periodicity parameter
    dg::Horner1d psit; //!< Toroidal flux \f$\psi_t (s)\f$
    dg::Horner1d psip; //!< Poloidal flux \f$\psi_p (s)\f$

    double R_0; //!< Major radius
    dg::RealFourier1d Raxis; //!< Magnetic axis \f$ R_a(\varphi)\f$ (same unit as \c R_0)
    dg::RealFourier1d Zaxis; //!< Magnetic axis \f$ Z_a(\varphi)\f$ (same unit as \c R_0)
    //dg::Horner2dRealFourier1d s; //!< Flux surface coordinate \f$ s(\bar R, \bar Z, \varphi)\f$ with \f$ \bar R = (R- R_a)/R_0\f$
    //dg::Horner2dRealFourier1d thetat; //!< Difference of flux angle coordinate to geometric poloidal angle coordinate \f$\tilde \theta_f \equiv (\theta + \lambda - \theta_g)(\bar R, \bar Z, \varphi)\f$
};


/*! @brief Convert a vmec file to Clebsch parameters (i.e. Cylindrical coordinates)
 *
 * using least squares optimisation
 * @tparam NcFile One of \c dg::file::NcFile classes. We use a template
 * parameter here so we do not directly depend on netcdf here.
 * @param vmec An open vmec file (or at least a file whose current group contains a copy of a vmec file)
 * @param eps Error tolerance for the least squares optimization
 */
template<class NcFile>
ClebschParameters vmec2feltor( const NcFile& vmec, double eps = 1e-8)
{
    ClebschParameters params;
    params.nfp = vmec.template get_var_as<int>( "nfp");
    params.R_0 = vmec.template get_var_as<double>( "Rmajor_p");
    unsigned n_tor = vmec.get_dim_size( "n_tor");
    params.Raxis = dg::RealFourier1d( vmec.template get_var_as<std::vector<double>>(
        "raxis_cc"), n_tor, 0, 2.*M_PI/(double)params.nfp); // 0 sine modes
    params.Zaxis = dg::RealFourier1d( vmec.template get_var_as<std::vector<double>>(
        "zaxis_cs", {1,n_tor-1}), 0, n_tor - 1, 2.*M_PI); // 0 cosine modes

    unsigned ns = vmec.template get_var_as<unsigned>("ns");
    auto phi = vmec.template get_var_as<std::vector<double>>( "phi"); // phi
    auto chi = vmec.template get_var_as<std::vector<double>>( "chi"); // chi

    auto chipf = vmec.template get_var_as<std::vector<double>>( "chipf"); // chi
    chipf[0] = 2.*chipf[1] - chipf[2];

    thrust::host_vector<double> s_grid(ns-1);
    for( unsigned u=0; u<ns; u++)
        s_grid[u] = (double)u/double(ns-1);

    double accuracy = 1e300;
    unsigned Ncoeff = 1;
    while( accuracy > eps && Ncoeff < ns)
    {
        std::vector<thrust::host_vector<double>> jac_poly1d( Ncoeff, s_grid);
        for( unsigned u=0; u<Ncoeff; u++)
            dg::blas1::transform( s_grid, jac_poly1d[u], [u](double s) { return pow(s,u+1);});

        auto coeffs = dg::least_squares( jac_poly1d, phi, 1.);
        coeffs.insert( coeffs.begin(), 0);

        dg::Horner1d horner( coeffs);
        thrust::host_vector<double> test( ns);
        dg::blas1::transform( s_grid, test, horner);
        dg::blas1::axpby( 1., phi, -1., test);
        accuracy = sqrt( dg::blas1::dot( test, test)/(double(ns-1)));
        params.psit = horner;
        std::cout << "Accuracy "<<accuracy<<"\n";
        std::cout << "Psi_t coeff "<<coeffs[0]<<" "<<coeffs[1]<<"\n";
        Ncoeff*=2;
    }
    double old_accuracy = 1e301;
    accuracy = 1e300;
    Ncoeff = 1;
    while( accuracy > eps && Ncoeff < ns && accuracy < old_accuracy)
    {
        old_accuracy = accuracy;
        std::vector<thrust::host_vector<double>> jac_poly1d( Ncoeff, s_grid);
        dg::blas1::transform( s_grid, jac_poly1d[0], [](double) { return 1.;});
        for( unsigned u=1; u<Ncoeff; u++)
        {
            //if( u==1)
            //    dg::blas1::evaluate( jac_poly1d[u], dg::equals(), [](double s){ return 1-s;}, s_grid);
            //else
            //    dg::blas1::evaluate( jac_poly1d[u], dg::equals(), [u](double s, double lm1, double lm2){
            //        double k = u-1;
            //        return ((2*k+1-s)*lm1 - k * lm2)/(k+1);
            //    }, s_grid, jac_poly1d[u-1], jac_poly1d[u-2]) ;
            dg::blas1::transform( s_grid, jac_poly1d[u], [u](double s) { return pow(s,u);});
        }

        auto coeffs = dg::least_squares( jac_poly1d, chi, 1.);

        //dg::Laguerre1d horner( coeffs);
        dg::Horner1d horner( coeffs);
        thrust::host_vector<double> test( ns);
        dg::blas1::transform( s_grid, test, horner);
        dg::blas1::axpby( 1., chi, -1., test);
        accuracy = sqrt( dg::blas1::dot( test, test)/(double(ns-1)));
        //params.psip = horner;
        std::cout << "N "<<Ncoeff<<"\n";
        std::cout << "Accuracy "<<accuracy<<"\n";
        std::cout << "psi_p(0) "<<horner(0)<<" psi_p (1) "<< horner(1)<<"\n";
        std::cout << "Psi_p coeff "<<coeffs[0]<<" "<<coeffs[1]<<" "<<coeffs[2]<<" "<<coeffs[3]<<"\n";
        Ncoeff+=1;
    }

    thrust::host_vector<unsigned> ns_grid( ns);
    for( unsigned u=0; u<ns_grid.size(); u++)
        ns_grid[u] = u;


    thrust::host_vector<double> theta_grid(100);
    for( unsigned u=0; u<theta_grid.size(); u++)
        theta_grid[u] = 2.*M_PI*u/double(theta_grid.size());

    thrust::host_vector<double> phi_grid(100/params.nfp);
    for( unsigned u=0; u<phi_grid.size(); u++)
        phi_grid[u] = 2.*M_PI*u/double(phi_grid.size()*params.nfp);

    unsigned mn_mode = vmec.get_dim_size( "mn_mode");
    std::vector<std::vector<double>> rmnc( ns), zmns( ns);
    std::vector<double> xm, xn;
    for( unsigned u=0; u<ns; u++)
    {
        vmec.get_var( "rmnc", {{u, 0}, {1, mn_mode}}, rmnc[u]);
        vmec.get_var( "zmns", {{u, 0}, {1, mn_mode}}, zmns[u]);
    }
    vmec.get_var( "xm", {0, mn_mode}, xm);
    vmec.get_var( "xn", {0, mn_mode}, xn);
    std::cout << "Kronecker\n";
    dg::Timer t;
    t.tic();
    auto rijk = dg::kronecker( [&](unsigned s, double t, double p){
            double rijk = 0;
            for( unsigned u=0; u<mn_mode; u++)
                rijk += rmnc[s][u]*cos( xm[u]*t - xn[u]*p);
            return rijk - params.Raxis(p); // we approximate in normalized coords
        }, ns_grid, theta_grid, phi_grid);
    t.toc();
    std::cout << "Took "<<t.diff()<<"\n";
    t.tic();
    auto zijk = dg::kronecker( [&](unsigned s, double t, double p){
            double zijk = 0;
            for( unsigned u=0; u<mn_mode; u++)
                zijk += zmns[s][u]*sin( xm[u]*t - xn[u]*p);
            return zijk - params.Zaxis(p);
        }, ns_grid, theta_grid, phi_grid);
    auto pijk = dg::kronecker( [&](unsigned s, double t, double p){
            return p;
        }, ns_grid, theta_grid, phi_grid);
    auto sijk = dg::kronecker( [&](unsigned s, double t, double p){
            return s_grid[s];
        }, ns_grid, theta_grid, phi_grid);
    t.toc();
    std::cout << "Took "<<t.diff()<<"\n";
    // least squares
    //
    t.tic();
    unsigned num_r = 8, num_z = 8, num_p = 8;
    std::vector<thrust::host_vector<double>> jac( num_r * num_z * 8, rijk);
    for( unsigned r=0; r<num_r; r++)
    for( unsigned z=0; z<num_z; z++)
    for( unsigned p=0; p<8; p++)
        dg::blas1::evaluate( jac[ (r*num_z + z)*8 +  p], dg::equals(),
         [&](double R, double Z, double P){
             if( p == 0)
                 return pow(R,r)*pow(Z,z);
             if( p < 5)
                 return pow(R,r)*pow(Z,z)*cos(p*params.nfp*P);
             return pow(R,r)*pow(Z,z)*sin((p-4)*params.nfp*P);
         }, rijk, zijk, pijk);
    t.toc();
    std::cout << "Jacobian assembly took "<<t.diff()<<"\n";

    auto coeffs = dg::least_squares( jac, sijk, 1.);
    dg::Horner2dRealFourier1d functor3d( coeffs, num_r, num_z, 5, 3, 2.*M_PI/params.nfp);
    auto test = rijk;
    dg::blas1::evaluate( test, dg::equals(), functor3d, rijk, zijk, pijk);
    dg::blas1::axpby( 1., sijk, -1., test);
    accuracy = sqrt( dg::blas1::dot( test, test)/(double(ns-1)));
    std::cout << "Accuracy "<<accuracy<<"\n";

    return params;

}
} //namespace vmec
} //namespace geo
} //namespace dg
