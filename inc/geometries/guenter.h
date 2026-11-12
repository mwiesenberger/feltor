#pragma once

#include <iostream>
#include <cmath>
#include <vector>

#include "dg/blas.h"

#include "magnetic_field.h"

/*!@file
 *
 * MagneticField objects
 */
namespace dg
{
namespace geo
{
/**
 * @brief Contains the Guenter type flux functions
 */
namespace guenter
{
///@addtogroup guenter
///@{


/**
 * @brief \f[\cos(\pi(R-R_0)/2)\cos(\pi Z/2)\f]
 */
struct Psip : public aCylindricalFunctor<Psip>
{
    /// @param R_0 major radius
    Psip(double R_0 ):   R_0(R_0) {}
    double do_compute(double R, double Z) const
    {
        return cos(M_PI*0.5*(R-R_0))*cos(M_PI*Z*0.5);
    }
  private:
    double R_0;
};
/**
 * @brief \f[-\pi\sin(\pi(R-R_0)/2)\cos(\pi Z/2)/2\f]
 */
struct PsipR : public aCylindricalFunctor<PsipR>
{
    /// @param R_0 major radius
    PsipR(double R_0 ):   R_0(R_0) {}
    double do_compute(double R, double Z) const
    {
        return -M_PI*0.5*sin(M_PI*0.5*(R-R_0))*cos(M_PI*Z*0.5);
    }
  private:
    double R_0;
};
/**
 * @brief \f[-\pi^2\cos(\pi(R-R_0)/2)\cos(\pi Z/2)/4\f]
 */
struct PsipRR : public aCylindricalFunctor<PsipRR>
{
    /// @param R_0 major radius
    PsipRR(double R_0 ):   R_0(R_0) {}
    double do_compute(double R, double Z) const
    {
        return -M_PI*M_PI*0.25*cos(M_PI*0.5*(R-R_0))*cos(M_PI*Z*0.5);
    }
  private:
    double R_0;
};
/**
 * @brief \f[-\pi\cos(\pi(R-R_0)/2)\sin(\pi Z/2)/2\f]
 */
struct PsipZ : public aCylindricalFunctor<PsipZ>

{
    /// @param R_0 major radius
    PsipZ(double R_0 ):   R_0(R_0) {}
    double do_compute(double R, double Z) const
    {
        return -M_PI*0.5*cos(M_PI*0.5*(R-R_0))*sin(M_PI*Z*0.5);
    }
  private:
    double R_0;
};
/**
 * @brief \f[-\pi^2\cos(\pi(R-R_0)/2)\cos(\pi Z/2)/4\f]
 */
struct PsipZZ : public aCylindricalFunctor<PsipZZ>
{
    /// @param R_0 major radius
    PsipZZ(double R_0 ):   R_0(R_0){}
    double do_compute(double R, double Z) const
    {
        return -M_PI*M_PI*0.25*cos(M_PI*0.5*(R-R_0))*cos(M_PI*Z*0.5);
    }
  private:
    double R_0;
};
/**
 * @brief \f[ \pi^2\sin(\pi(R-R_0)/2)\sin(\pi Z/2)/4\f]
 */
struct PsipRZ : public aCylindricalFunctor<PsipRZ>
{
    /// @param R_0 major radius
    PsipRZ(double R_0 ):   R_0(R_0) {}
    double do_compute(double R, double Z) const
    {
        return M_PI*M_PI*0.25*sin(M_PI*0.5*(R-R_0))*sin(M_PI*Z*0.5);
    }
  private:
    double R_0;
};

/**
 * @brief \f[ I_0\f]
 */
struct Ipol : public aCylindricalFunctor<Ipol>
{
    /// @param I_0 the constant current
    Ipol( double I_0):   I_0(I_0) {}
    double do_compute(double, double) const { return I_0; }
    private:
    double I_0;
};
/**
 * @brief \f[0\f]
 */
struct IpolR : public aCylindricalFunctor<IpolR>
{
    ///@brief No parameters
    IpolR(  ) {}
    double do_compute(double, double) const { return 0; }
    private:
};
/**
 * @brief \f[0\f]
 */
struct IpolZ : public aCylindricalFunctor<IpolZ>
{
    ///@brief No parameters
    IpolZ(  ) {}
    double do_compute(double, double) const { return 0; }
    private:
};

/**
 * @brief \f$ \psi_p\f$ and its first and second derivatives
 * @param R_0 major radius
 * @return \c Psip, \c PsipR, \c PsipZ, \c PsipRR, \c PsipRZ, \c PsipZZ
 */
inline CylindricalFunctorsLvl2 createPsip( double R_0)
{
    return CylindricalFunctorsLvl2( Psip(R_0), PsipR(R_0), PsipZ(R_0),
        PsipRR(R_0), PsipRZ(R_0), PsipZZ(R_0));
}
/**
 * @brief \f$ I\f$ and its first derivatives
 * @param I_0 the constant current
 * @return \c Ipol, \c IpolR, \c IpolZ
 */
inline CylindricalFunctorsLvl1 createIpol( double I_0)
{
    return CylindricalFunctorsLvl1( Ipol(I_0), IpolR(), IpolZ());
}
///@}
} //namespace guenter

/**
 * @brief Create a Guenter Magnetic field

 * \f[\psi_p(R,Z) = \cos(\pi(R-R_0)/2)\cos(\pi Z/2),\quad I(\psi_p) = I_0\f]
 * @param R_0 the major radius
 * @param I_0 the current
 * @return A magnetic field object
 * @ingroup guenter
 */
inline dg::geo::TokamakMagneticField createGuenterField( double R_0, double I_0)
{
    MagneticFieldParameters params = { 1., 1., 0.,
            equilibrium::guenter, modifier::none, description::square};
    return TokamakMagneticField( R_0, guenter::createPsip(R_0), guenter::createIpol(I_0), params);
}
} //namespace geo
}//namespace dg
