#pragma once

#include <string>
#include <vector>
#include <functional>

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"

#include "thermal.h"
#include "parameters.h"

#include "init.h"
#include "../feltor/feltordiag.h" // for static lists
#include "../feltor/common.h"

namespace thermal{

// This file constitutes the diagnostics module for thermal
// The way it works is that it generates global lists of Records that describe what goes into the file
// You can register you own diagnostics in one of three diagnostics lists (static 3d, dynamic 3d and
// dynamic 2d) further down
// which will then be applied during a simulation
struct BPerp{
    //b_perp
    DG_DEVICE void operator()(double A,
        double d0A, double d1A,
        double& bp0, double& bp1, //bperp
        double b_2,
        double curvKappa0,  double curvKappa1
        ){
        bp0 = (b_2*d1A  + A*curvKappa0);
        bp1 = (- b_2*d0A + A*curvKappa1);
    }
};


//From here on, we use the typedefs to ease the notation

struct Variables{
    thermal::Explicit<dg::x::CylindricalGrid3d, dg::x::IDMatrix, dg::x::DMatrix, dg::x::DVec>& f;
    std::array<std::vector<dg::x::DVec>,6>& y0;
    thermal::Parameters p;
    dg::geo::TokamakMagneticField mag;
    const std::array<dg::x::DVec, 2>& gradPsip;
    std::array<dg::x::DVec, 2> tmp;
    std::array<dg::x::DVec, 2> tmp2;
    std::array<dg::x::DVec, 2> tmp3;
    double duration;
    const unsigned* nfailed;
};

struct Record{
    bool species_dependent; // whether variable should be prepended with species name
    std::string name;
    std::string long_name;
    bool integral; //indicates whether the function should be time-integrated
    std::function<void( dg::x::DVec&, Variables&, unsigned)> function;
};

///%%%%%%%%%%%%%%%%%%%%EXTEND LISTS WITH YOUR DIAGNOSTICS HERE%%%%%%%%%%%%%%%%%%%%%%
///%%%%%%%%%%%%%%%%%%%%EXTEND LISTS WITH YOUR DIAGNOSTICS HERE%%%%%%%%%%%%%%%%%%%%%%
///%%%%%%%%%%%%%%%%%%%%EXTEND LISTS WITH YOUR DIAGNOSTICS HERE%%%%%%%%%%%%%%%%%%%%%%
// Here is a list of useful 1d variables of general interest
std::vector<dg::file::Record<double(Variables&), dg::file::LongNameAttribute>> diagnostics1d_list = {
    {"failed", "Accumulated Number of failed steps",
        []( Variables& v ) {
            return *v.nfailed;
        }
    },
    {"duration", "Computation time between the latest 3d outputs (without the output time itself)",
        []( Variables& v ) {
            return v.duration;
        }
    },
    {"nsteps", "Accumulated Number of calls to the right-hand-side (including failed steps)",
        [](Variables& v ) {
            return v.f.called();
        }
    }
};
//Here is a list of static (time-independent) 3d variables that go into the output
//Cannot be thermal internal variables
std::vector<dg::file::Record<void( dg::x::HVec&, const dg::geo::TokamakMagneticField&, const dg::x::CylindricalGrid3d&), dg::file::LongNameAttribute>> diagnostics3d_static_list =
    feltor::diagnostics3d_static_list;

// Here are all 3d outputs we want to have
std::vector<Record> diagnostics3d_list = { // 2 + 6*s
    {true, "n", "gyro-centre density", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("N", s), result);
        }
    },
    {true, "tperp", "perpendicular temperature", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("Tperp", s), result);
        }
    },
    {true, "tpara", "parallel temperature", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("Tpara", s), result);
        }
    },
    {true, "u", "parallel velocity", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("U", s), result);
        }
    },
    {true, "uperp", "perpendicular heat flux velocity", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("Uperp", s), result);
        }
    },
    {true, "upara", "parallel heat flux velocity", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("Upara", s), result);
        }
    },
    {false, "phi", "electric potential", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::copy(v.f.potential(), result);
        }
    },
    {false, "aparallel", "parallel magnetic potential", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::copy(v.f.apar(), result);
        }
    }
};

std::vector<dg::file::Record<void(dg::x::HVec&, Variables&, const dg::x::CylindricalGrid3d&), dg::file::LongNameAttribute>> diagnostics2d_static_list = {
    { "Psip2d", "Flux-function psi",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.psip(), grid);
        }
    },
    { "PsipR2d", "Flux-function psi R-derivative",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.psipR(), grid);
        }
    },
    { "PsipZ2d", "Flux-function psi Z-derivative",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.psipZ(), grid);
        }
    },
    { "PsipRR2d", "Flux-function psi RR-derivative",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.psipRR(), grid);
        }
    },
    { "PsipRZ2d", "Flux-function psi RZ-derivative",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.psipRZ(), grid);
        }
    },
    { "PsipZZ2d", "Flux-function psi ZZ-derivative",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.psipZZ(), grid);
        }
    },
    { "Ipol", "Poloidal current",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.ipol(), grid);
        }
    },
    { "IpolR", "Poloidal current R-derivative",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.ipolR(), grid);
        }
    },
    { "IpolZ", "Poloidal current Z-derivative",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( v.mag.ipolZ(), grid);
        }
    },
    { "Rho_p", "Normalized Poloidal flux label",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( dg::geo::RhoP( v.mag), grid);
        }
    },
    { "Bmodule", "Magnetic field strength",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            result = dg::pullback( dg::geo::Bmodule(v.mag), grid);
        }
    },
    { "Divb", "The divergence of the magnetic unit vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            dg::assign(  dg::pullback(dg::geo::Divb(v.mag), grid), result);
        }
    },
    { "InvB", "Inverse of Bmodule",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid ){
            dg::assign(  dg::pullback(dg::geo::InvB(v.mag), grid), result);
        }
    },
    { "CurvatureKappaR", "R-component of the Kappa B curvature vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& ){
            dg::assign( v.f.perp().curvKappa()[0], result);
        }
    },
    { "CurvatureKappaZ", "Z-component of the Kappa B curvature vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& ){
            dg::assign( v.f.perp().curvKappa()[1], result);
        }
    },
    { "CurvatureKappaP", "Contravariant Phi-component of the Kappa B curvature vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& ){
            dg::assign( v.f.perp().curvKappa()[2], result);
        }
    },
    { "DivCurvatureKappa", "Divergence of the Kappa B curvature vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& ){
            dg::assign( v.f.perp().divCurvKappa(), result);
        }
    },
    { "CurvatureNablaR", "R-component of the Nabla B curvature vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& ){
            dg::assign( v.f.perp().curvNabla()[0], result);
        }
    },
    { "CurvatureNablaZ", "Z-component of the Nabla B curvature vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& ){
            dg::assign( v.f.perp().curvNabla()[1], result);
        }
    },
    { "bphi", "Covariant Phi-component of the magnetic unit vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& ){
            dg::assign( v.f.perp().bphi(), result);
        }
    },
    { "BHatR", "R-component of the magnetic field unit vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid){
            result = dg::pullback( dg::geo::BHatR(v.mag), grid);
        }
    },
    { "BHatZ", "Z-component of the magnetic field unit vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid){
            result = dg::pullback( dg::geo::BHatZ(v.mag), grid);
        }
    },
    { "BHatP", "P-component of the magnetic field unit vector",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid){
            result = dg::pullback( dg::geo::BHatP(v.mag), grid);
        }
    },
    { "NormGradPsip", "Norm of gradient of Psip",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& grid){
            result = dg::pullback(
                dg::geo::SquareNorm( dg::geo::createGradPsip(v.mag),
                    dg::geo::createGradPsip(v.mag)), grid);
        }
    },
    { "Wall", "Wall Region",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d&  ){
            dg::assign( v.f.sources().get_wall(), result);
        }
    },
    { "Sheath", "Sheath Region",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d&  ){
            dg::assign( v.f.para().get_sheath(), result);
        }
    },
    { "SheathCoordinate", "Sheath Coordinate of field lines",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d&  ){
            dg::assign( v.f.para().get_sheath_coordinate(), result);
        }
    },
    { "vol2d", "Volume form (including R) in 2d",
        []( dg::x::HVec& result, Variables&, const dg::x::CylindricalGrid3d& grid ){
            result = dg::create::volume(grid);
        }
    }
};

std::vector<dg::file::Record<void(dg::x::HVec&, Variables&, const dg::x::CylindricalGrid3d&, unsigned s), dg::file::LongNameAttribute>> diagnostics2d_static_init = {
    { "Nprof", "Density profile (that the source may force)",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.f.sources().get_source_prof(0,s), result);
        }
    },
    { "Pperpprof", "Pperp profile (that the source may force)",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.f.sources().get_source_prof(1,s), result);
        }
    },
    { "Pparaprof", "Ppara profile (that the source may force)",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.f.sources().get_source_prof(2,s), result);
        }
    },
    { "SN", "Density source profile (influx)",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.f.sources().get_source_region(0,s), result);
        }
    },
    { "SPperp", "Pperp source profile (influx)",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.f.sources().get_source_region(1,s), result);
        }
    },
    { "SPpara", "Ppara source profile (influx)",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.f.sources().get_source_region(2,s), result);
        }
    },
    { "Ninit", "Initial density condition",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.y0[0][s], result);
        }
    },
    { "Pperpinit", "Initial perp pressure condition",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.y0[1][s], result);
        }
    },
    { "Pparainit", "Initial para pressure condition",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.y0[2][s], result);
        }
    },
    { "Winit", "Initial canonical velocity condition",
        []( dg::x::HVec& result, Variables& v, const dg::x::CylindricalGrid3d& , unsigned s ){
            dg::assign( v.y0[3][s], result);
        }
    },
};

// and here are all the 2d outputs we want to produce (currently ~ 150)
// Call within species loop after updateQuantities
std::vector<Record> basicDiagnostics2d_list = { // 22
    {true, "n", "gyro-centre density", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("N",s), result);
        }
    },
    {true, "tperp", "perpendicular temperature", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("Tperp",s), result);
        }
    },
    {true, "tpara", "parallel temperature", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("Tpara",s), result);
        }
    },
    {true, "u", "parallel velocity", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("U",s), result);
            // U = W - q/m Apar
        }
    },
    {true, "uperp", "perpendicular heat flux velocity Qperp/Pperp", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("Uperp",s), result);
        }
    },
    {true, "upara", "parallel heat flux velocity Qpara/Ppara", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            dg::blas1::copy(v.f.get("Upara",s), result);
        }
    },
    {false, "phi", "electric potential", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::copy(v.f.potential(), result);
        }
    },
    {false, "aparallel", "parallel magnetic potential", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::copy(v.f.apar(), result);
        }
    },
    {true, "psi0", "Potential 0", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
             dg::blas1::copy(v.f.get("Psi0",s), result);
        }
    },
    {true, "psi1", "Potential 0", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
             dg::blas1::copy(v.f.get("Psi1",s), result);
        }
    },
    {true, "psi2", "Potential 0", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
             dg::blas1::copy(v.f.get("Psi2",s), result);
        }
    },
    {true, "psi3", "Potential 0", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
             dg::blas1::copy(v.f.get("Psi3",s), result);
        }
    },
    {true, "gammaN", "Adjoint Gamma N", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
             dg::blas1::copy(v.f.solvers().gammaN(s), result);
        }
    },
    /// -----------------Miscellaneous additions --------------------//
    {false, "vorticity", "Minus Lap_perp of potential", false,
        []( dg::x::DVec& result, Variables& v , unsigned) {
            // has no jump terms
            v.f.solvers().compute_lapMperpP(v.f.potential(), result);
        }
    },
    {true, "laplace_n", "Positive Lap_perp of density", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            v.f.perp().compute_perp_laplace(-1.0, v.f.get("N", s), 1, v.tmp[0], v.tmp[1], 0., result);
        }
    },
    // Does not work due to direct application of Laplace
    // The Laplacian of Aparallel looks smooth in paraview
    {false, "apar_vorticity", "Minus Lap_perp of magnetic potential", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            v.f.perp().compute_perp_laplace(-1.0, v.f.apar(), 1, v.tmp[0], v.tmp[1], 0., result);
        }
    },
    {true, "dssu", "2nd parallel derivative of velocity", false,
        []( dg::x::DVec& result, Variables& v, unsigned s ) {
            v.f.para().fieldaligned()( dg::geo::einsMinus, v.f.get("U", s), v.tmp[0]);
            v.f.para().fieldaligned()( dg::geo::zeroForw,  v.f.get("U", s), v.tmp2[0]);
            v.f.para().fieldaligned()( dg::geo::einsPlus,  v.f.get("U", s), v.tmp3[0]);
            v.f.para().update_parallel_bc_2nd( v.f.para().fieldaligned(),
                v.tmp[0], v.tmp[1], v.tmp[2], v.p.bcx, 0.);
            dg::geo::dss_centered( v.f.para().fieldaligned(), 1.,
                v.tmp[0], v.tmp[1], v.tmp[2], 0, result);

        }
    },
    {false, "lperpinv", "Perpendicular (electron) density gradient length scale", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::pointwiseDot( 1., v.f.get( "dxF N", 0), v.f.get( "dxF N", 0),
                                     1., v.f.get( "dyF N", 0), v.f.get( "dyF N", 0),
                                     0., result);
            dg::blas1::pointwiseDivide( result, v.f.get( "N", 0), result);
            dg::blas1::pointwiseDivide( result, v.f.get( "N", 0), result);
            // ((grad N)/N)**2
            dg::blas1::transform( result, result, dg::SQRT<double>());
        }
    },
    {false, "perpaligned", "Perpendicular (electron) density alignement", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::pointwiseDot( 1., v.f.get( "dxF N", 0), v.f.get( "dxF N", 0),
                                     1., v.f.get( "dyF N", 0), v.f.get( "dyF N", 0),
                                     0., result);
            dg::blas1::pointwiseDivide( result, v.f.get( "N", 0), result);
        }
    },
    {false, "lparallelinv", "Parallel (electron) density gradient length scale", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::geo::ds_centered( v.f.para().fieldaligned(), 1.,
                v.f.get( "N -1", 0), v.f.get( "N +1", 0), 0, result);
            dg::blas1::pointwiseDivide( result, v.f.get("N", 0), result);
            dg::blas1::pointwiseDot ( result, result, result);
            dg::blas1::transform( result, result, dg::SQRT<double>());
        }
    },
    {false, "aligned", "Parallel (electron) density alignement", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::geo::ds_centered( v.f.para().fieldaligned(), 1.,
                v.f.get( "N -1", 0), v.f.get( "N +1", 0), 0, result);
            dg::blas1::pointwiseDot ( result, result, result);
            dg::blas1::pointwiseDivide( result, v.f.get("N", 0), result);
        }
    },
    /// ------------------ Correlation terms --------------------//
    {false, "ne2", "Square of electron density", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::pointwiseDot(
                v.f.get("N", 0), v.f.get("N", 0), result);
        }
    },
    {false, "phi2", "Square of electron potential", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::pointwiseDot(
                v.f.potential(), v.f.potential(), result);
        }
    },
    {false, "nephi", "Product of electron potential and electron density", false,
        []( dg::x::DVec& result, Variables& v, unsigned ) {
            dg::blas1::pointwiseDot(
                v.f.potential(), v.f.get("N", 0), result);
        }
    }
};

// TODO Fix conservation theorems and Probes
/*
std::vector<Record> MassConsDiagnostics2d_list = { // 26
    /// ------------------ Density terms ------------------------//
    ////////////////// electron particle flux /////////////////////
    {true, "jsnE_tt", "Radial particle flux: ExB contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            // ExB Dot GradPsi
            routines::jacobian( v.f.bhatgB(), v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( result, v.f.density(0), result);
        }
    },
    {true, "divnE_tt", "Divergence of ExB particle flux (Time average)", true,
        []( dg::x::DVec& result, Variables& v) {
            routines::dot( v.f.curv(), v.f.gradP(0), result);
            dg::blas1::pointwiseDot( result, v.f.density(0), result);
            routines::jacobian( 1., v.f.bhatgB(), v.f.gradP(0), v.f.gradN(0), 1., result);
        }
    },
    {true, "jscurvn_tt", "Radial particle flux: curvature contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            routines::dot( v.f.curv(), v.gradPsip, result);
            dg::blas1::pointwiseDot( v.p.tau[0], v.f.density(0), result, 0., result);
        }
    },
    {true, "divcurvn_tt", "Divergence of curvature term (Time average)", true,
        []( dg::x::DVec& result, Variables& v) {
            routines::dot( v.p.tau[0], v.f.curv(), v.f.gradN(0), 0., result);
        }
    },
    {true, "jscurvkappan_tt", "Radial particle flux: curvature contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot( v.p.mu[0], v.f.density(0), v.f.velocity(0),
                    v.f.velocity(0), 0., result);
            routines::dot( v.f.curvKappa(), v.gradPsip, v.tmp[0]);
            dg::blas1::pointwiseDot( v.tmp[0], result, result);
        }
    },
    {true, "divcurvkappan_tt", "Divergence of curvature term (Time average)", true,
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::pointwiseDot( v.p.mu[0], v.f.density(0), v.f.velocity(0),
                    v.f.velocity(0), 0., v.tmp3[0]);
            v.f.centered_div( v.tmp3[0], v.f.curvKappa(), v.tmp2[0], result);
        }
    },
    {true, "jsnA_tt", "Radial particle flux: magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_bperp(v.tmp);
            routines::dot( v.tmp, v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., v.f.density(0), v.f.velocity(0), result, 0., result);
        }
    },
    {true, "divnA_tt", "Divergence of magnetic flutter particle flux (Time average)", true,
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::pointwiseDot( v.f.density(0), v.f.velocity(0), v.tmp3[0]);
            v.f.compute_bperp(v.tmp);
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {true, "jsdia_tt", "Radial particle flux: diamagnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            // u_D Dot GradPsi
            routines::jacobian( v.p.tau[0], v.f.bhatgB(), v.f.gradN(0), v.gradPsip, 0., result);
        }
    },
    {true, "lnperp_tt", "Perpendicular diffusion (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_perp_diffusiveN( 1., v.f.density(0), v.tmp[0],
                    v.tmp[1], 0., result);
        }
    },
    {true, "lnparallel_tt", "Parallel diffusion (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::axpby( v.p.nu_parallel_n, v.f.lapParN(0), 0., result);
        }
    },
    {true, "sn_tt", "Source term for density (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::copy( v.f.density_source(0), result);
        }
    },
    {true, "divjnpar_tt", "Divergence of Parallel velocity term for density (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::copy( v.f.divNUb(0), result);
        }
    }
};

std::vector<Record> EnergyDiagnostics2d_list = { // 23
    /// ------------------- Energy terms ------------------------//
    {true, "nlnn", "Entropy", false,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::transform( v.f.density(0), result, feltor::PositiveLN());
            dg::blas1::pointwiseDot( result, v.f.density(0), result);
        }
    },
    {false, "aperp2", "Magnetic energy", false,
        []( dg::x::DVec& result, Variables& v ) {
            if( v.p.beta == 0) // avoid divison by zero
            {
                dg::blas1::scal( result, 0.);
            }
            else
            {
                dg::tensor::scalar_product3d( 1./2./v.p.beta, 1.,
                    v.f.gradA()[0], v.f.gradA()[1], v.f.gradA()[2],
                    v.f.projection(), 1., //grad_perp
                    v.f.gradA()[0], v.f.gradA()[1], v.f.gradA()[2], 0., result);
            }
        }
    },
    {false, "ue2", "ExB energy", false,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot( 0.5, v.f.density(1), v.f.uE2(), 0., result);
        }
    },
    {true, "nu2", "Parallel kinetic energy", false,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot( -0.5*v.p.mu[0], v.f.density(0),
                v.f.velocity(0), v.f.velocity(0), 0., result);
        }
    },
    /// ------------------- Energy dissipation ----------------------//
    {"resistivity_tt", "Energy dissipation through resistivity (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot( 1., v.f.velocity(1), v.f.density(0), -1.,
                    v.f.velocity(0), v.f.density(0), 0., v.tmp[0]);
            dg::blas1::pointwiseDot( 1., v.f.velocity(1), v.f.density(1), -1.,
                    v.f.velocity(0), v.f.density(0), 0., v.tmp[1]);
            dg::blas1::pointwiseDot( -v.p.eta, v.tmp[0], v.tmp[1], 0., result);
        }
    },
    {"see_tt", "Energy sink/source for electrons", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::evaluate( result, dg::equals(),
                RadialEnergyFlux( v.p.tau[0], v.p.mu[0], -1.),
                v.f.density(0), v.f.velocity(0), v.f.potential(0),
                v.f.density_source(0)
            );
        }
    },
    {"sei_tt", "Energy sink/source for ions", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::evaluate( result, dg::equals(),
                RadialEnergyFlux( v.p.tau[1], v.p.mu[1], 1.),
                v.f.density(1), v.f.velocity(1), v.f.potential(1),
                v.f.density_source(1)
            );
        }
    },
    /// ------------------ Energy flux terms ------------------------//
    {"jsee_tt", "Radial electron energy flux without magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::subroutine(
                RadialEnergyFlux( v.p.tau[0], v.p.mu[0], -1.),
                v.f.density(0), v.f.velocity(0), v.f.potential(0),
                v.f.gradP(0)[0], v.f.gradP(0)[1], v.f.gradP(0)[2],
                v.tmp[0], v.tmp[1], v.tmp[2],
                v.f.bhatgB()[0], v.f.bhatgB()[1], v.f.bhatgB()[2],
                v.f.curv()[0], v.f.curv()[1], v.f.curv()[2],
                v.f.curvKappa()[0], v.f.curvKappa()[1], v.f.curvKappa()[2]
            );
            routines::dot( v.tmp, v.gradPsip, result);
        }
    },
    {"jseea_tt", "Radial electron energy flux: magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::subroutine(
                RadialEnergyFlux( v.p.tau[0], v.p.mu[0], -1.),
                v.f.density(0), v.f.velocity(0), v.f.potential(0), v.f.apar(),
                v.f.gradA()[0], v.f.gradA()[1], v.f.gradA()[2],
                v.tmp[0], v.tmp[1], v.tmp[2],
                v.f.bhatgB()[0], v.f.bhatgB()[1], v.f.bhatgB()[2],
                v.f.curvKappa()[0], v.f.curvKappa()[1], v.f.curvKappa()[2]
            );
            routines::dot( v.tmp, v.gradPsip, result);
        }
    },
    {"jsei_tt", "Radial ion energy flux without magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::subroutine(
                RadialEnergyFlux( v.p.tau[1], v.p.mu[1], 1.),
                v.f.density(1), v.f.velocity(1), v.f.potential(1),
                v.f.gradP(1)[0], v.f.gradP(1)[1], v.f.gradP(1)[2],
                v.tmp[0], v.tmp[1], v.tmp[2],
                v.f.bhatgB()[0], v.f.bhatgB()[1], v.f.bhatgB()[2],
                v.f.curv()[0], v.f.curv()[1], v.f.curv()[2],
                v.f.curvKappa()[0], v.f.curvKappa()[1], v.f.curvKappa()[2]
            );
            routines::dot( v.tmp, v.gradPsip, result);
        }
    },
    {"jseia_tt", "Radial ion energy flux: magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::subroutine(
                RadialEnergyFlux( v.p.tau[1], v.p.mu[1], 1.),
                v.f.density(1), v.f.velocity(1), v.f.potential(1), v.f.apar(),
                v.f.gradA()[0], v.f.gradA()[1], v.f.gradA()[2],
                v.tmp[0], v.tmp[1], v.tmp[2],
                v.f.bhatgB()[0], v.f.bhatgB()[1], v.f.bhatgB()[2],
                v.f.curvKappa()[0], v.f.curvKappa()[1], v.f.curvKappa()[2]
            );
            routines::dot( v.tmp, v.gradPsip, result);
        }
    },
    {"divee_tt", "Radial electron energy flux without magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::subroutine(
                RadialEnergyFlux( v.p.tau[0], v.p.mu[0], -1.),
                v.f.density(0), v.f.velocity(0), v.f.potential(0),
                v.f.gradP(0)[0], v.f.gradP(0)[1], v.f.gradP(0)[2],
                v.tmp[0], v.tmp[1], v.tmp[2],
                v.f.bhatgB()[0], v.f.bhatgB()[1], v.f.bhatgB()[2],
                v.f.curv()[0], v.f.curv()[1], v.f.curv()[2],
                v.f.curvKappa()[0], v.f.curvKappa()[1], v.f.curvKappa()[2]
            );
            v.f.centered_div( 1., v.tmp, v.tmp2[0], result);
        }
    },
    {"diveea_tt", "Radial electron energy flux: magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::subroutine(
                RadialEnergyFlux( v.p.tau[0], v.p.mu[0], -1.),
                v.f.density(0), v.f.velocity(0), v.f.potential(0), v.f.apar(),
                v.f.gradA()[0], v.f.gradA()[1], v.f.gradA()[2],
                v.tmp[0], v.tmp[1], v.tmp[2],
                v.f.bhatgB()[0], v.f.bhatgB()[1], v.f.bhatgB()[2],
                v.f.curvKappa()[0], v.f.curvKappa()[1], v.f.curvKappa()[2]
            );
            v.f.centered_div( 1., v.tmp, v.tmp2[0], result);
        }
    },
    {"divei_tt", "Radial ion energy flux without magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::subroutine(
                RadialEnergyFlux( v.p.tau[1], v.p.mu[1], 1.),
                v.f.density(1), v.f.velocity(1), v.f.potential(1),
                v.f.gradP(1)[0], v.f.gradP(1)[1], v.f.gradP(1)[2],
                v.tmp[0], v.tmp[1], v.tmp[2],
                v.f.bhatgB()[0], v.f.bhatgB()[1], v.f.bhatgB()[2],
                v.f.curv()[0], v.f.curv()[1], v.f.curv()[2],
                v.f.curvKappa()[0], v.f.curvKappa()[1], v.f.curvKappa()[2]
            );
            v.f.centered_div( 1., v.tmp, v.tmp2[0], result);
        }
    },
    {"diveia_tt", "Radial ion energy flux: magnetic contribution (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::subroutine(
                RadialEnergyFlux( v.p.tau[1], v.p.mu[1], 1.),
                v.f.density(1), v.f.velocity(1), v.f.potential(1), v.f.apar(),
                v.f.gradA()[0], v.f.gradA()[1], v.f.gradA()[2],
                v.tmp[0], v.tmp[1], v.tmp[2],
                v.f.bhatgB()[0], v.f.bhatgB()[1], v.f.bhatgB()[2],
                v.f.curvKappa()[0], v.f.curvKappa()[1], v.f.curvKappa()[2]
            );
            v.f.centered_div( 1., v.tmp, v.tmp2[0], result);
        }
    },
    /// ------------------------ Energy dissipation terms ------------------//
    {"leeperp_tt", "Perpendicular electron energy dissipation (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_perp_diffusiveN( 1., v.f.density(0), v.tmp3[0], v.tmp3[1],
                    0., v.tmp[0]);
            v.f.compute_perp_diffusiveU( 1., v.f.velocity(0), v.f.density(0),
                    v.tmp3[0], v.tmp3[1], v.tmp3[2], v.tmp2[0], 0., v.tmp[1]);
            dg::blas1::evaluate( result, dg::equals(),
                RadialEnergyFlux( v.p.tau[0], v.p.mu[0], -1.),
                v.f.density(0), v.f.velocity(0), v.f.potential(0),
                v.tmp[0], v.tmp[1]
            );
        }
    },
    {"leiperp_tt", "Perpendicular ion energy dissipation (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_perp_diffusiveN( 1., v.f.density(1), v.tmp3[0], v.tmp3[1],
                    0., v.tmp[0]);
            v.f.compute_perp_diffusiveU( 1., v.f.velocity(1), v.f.density(1),
                    v.tmp3[0], v.tmp3[1], v.tmp3[2], v.tmp2[0], 0., v.tmp[1]);
            dg::blas1::evaluate( result, dg::equals(),
                RadialEnergyFlux( v.p.tau[1], v.p.mu[1], 1.),
                v.f.density(1), v.f.velocity(1), v.f.potential(1),
                v.tmp[0], v.tmp[1]
            );
        }
    },
    {"leeparallel_tt", "Parallel electron energy dissipation (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_parallel_diffusiveN( 0, v.tmp[0]);
            v.f.compute_parallel_diffusiveU( 0, v.tmp[1]);
            dg::blas1::evaluate( result, dg::equals(),
                RadialEnergyFlux( v.p.tau[0], v.p.mu[0], -1.),
                v.f.density(0), v.f.velocity(0), v.f.potential(0),
                v.tmp[0], v.tmp[1]
            );
        }
    },
    {"leiparallel_tt", "Parallel ion energy dissipation (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_parallel_diffusiveN( 1, v.tmp[0]);
            v.f.compute_parallel_diffusiveU( 1, v.tmp[1]);
            dg::blas1::evaluate( result, dg::equals(),
                RadialEnergyFlux( v.p.tau[1], v.p.mu[1], 1.),
                v.f.density(1), v.f.velocity(1), v.f.potential(1),
                v.tmp[0], v.tmp[1]
            );
        }
    },
    {"divjeepar_tt", "Divergence of Parallel electron energy flux (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            // Multiply out divNUb to get implementable form
            double z = -1.;
            dg::blas1::evaluate( result, dg::equals(),
                RadialEnergyFlux( v.p.tau[0], v.p.mu[0], z),
                v.f.density(0), v.f.velocity(0), v.f.potential(0),
                v.f.divNUb(0), 0.
            );
            dg::blas1::pointwiseDot( z*v.p.tau[0], v.f.velocity(0), v.f.dsN(0),
                    1., result);
            dg::blas1::pointwiseDot( z, v.f.velocity(0), v.f.density(0),
                    v.f.dsP(0), 1., result);
            dg::blas1::pointwiseDot( 1., v.f.density(0), v.f.velocity(0),
                    v.f.velocity(0), 0., v.tmp[0]);
            dg::blas1::pointwiseDot( z*v.p.mu[0], v.tmp[0], v.f.dsU(0), 1., result);
        }
    },
    {"divjeipar_tt", "Divergence of Parallel ion energy flux (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            // Multiply out divNUb to get implementable form
            double z = +1.;
            dg::blas1::evaluate( result, dg::equals(),
                RadialEnergyFlux( v.p.tau[1], v.p.mu[1], z),
                v.f.density(1), v.f.velocity(1), v.f.potential(1),
                v.f.divNUb(1), 0.
            );
            dg::blas1::pointwiseDot( z*v.p.tau[1], v.f.velocity(1), v.f.dsN(1),
                    1., result);
            dg::blas1::pointwiseDot( z, v.f.velocity(1), v.f.density(1),
                    v.f.dsP(1), 1., result);
            dg::blas1::pointwiseDot( 1., v.f.density(1), v.f.velocity(1),
                    v.f.velocity(1), 0., v.tmp[1]);
            dg::blas1::pointwiseDot( z*v.p.mu[1], v.tmp[1], v.f.dsU(1), 1., result);
        }
    }
};

std::vector<Record> ToroidalExBDiagnostics2d_list = { //27
    /// ------------------------ Vorticity terms ---------------------------//
    /// ----------------------with ion density -------------------------///
    {"oexbi", "ExB vorticity term with ion density", false,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDot( v.p.mu[1], result, v.f.density(1), 0., result);
        }
    },
    {"odiai", "Diamagnetic vorticity term with ion density", false,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( v.f.gradN(1), v.gradPsip, result);
            dg::blas1::scal( result, v.p.mu[1]*v.p.tau[1]);
        }
    },
    {"jsoexbi_tt", "ExB vorticity flux term with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB Dot GradPsi
            routines::jacobian( v.f.bhatgB(), v.f.gradP(0), v.gradPsip, v.tmp[0]);

            // Omega_E
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDot( v.p.mu[1], result, v.f.density(1), 0., result);

            // Multiply everything
            dg::blas1::pointwiseDot( 1., result, v.tmp[0], 0., result);
        }
    },
    {"divoexbi_tt", "ExB vorticity flux term with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB Dot GradPsi
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp); //u_E

            // Omega_E
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDot( v.p.mu[1], result, v.f.density(1), 0., v.tmp3[0]);

            // Divergence
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"jsoexbiUD_tt", "ExB vorticity flux term by diamagnetic velocity with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // bxGradN/B Dot GradPsi
            routines::jacobian( v.p.tau[1], v.f.bhatgB(), v.f.gradN(1), v.gradPsip, 0., result);

            // m Omega_E,phi
            routines::dot( v.f.gradP(0), v.gradPsip, v.tmp[0]);
            dg::blas1::pointwiseDot( v.p.mu[1], v.tmp[0], v.f.binv(), v.f.binv(), 0., v.tmp[0]);

            // Multiply everything
            dg::blas1::pointwiseDot( 1., result, v.tmp[0], 0., result);
        }
    },
    {"divoexbiUD_tt", "ExB vorticity flux term by diamagnetic velocity with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            routines::times( v.f.bhatgB(), v.f.gradN(1), v.tmp);

            // Omega_E
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( v.p.tau[1]*v.p.mu[1], result, v.f.binv(), v.f.binv(), 0., v.tmp3[0]);

            // Divergence
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"jsodiaiUE_tt", "Diamagnetic vorticity flux by ExB veloctiy term with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB Dot GradPsi
            routines::jacobian( v.f.bhatgB(), v.f.gradP(0), v.gradPsip, result);

            // Omega_D,phi
            routines::dot( v.p.mu[1]*v.p.tau[1], v.f.gradN(1), v.gradPsip, 0., v.tmp[0]);

            // Multiply everything
            dg::blas1::pointwiseDot( 1., result, v.tmp[0], 0., result);
        }
    },
    {"divodiaiUE_tt", "Diamagnetic vorticity flux by ExB veloctiy term with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB
            routines::times( v.f.bhatgB(), v.f.gradP(0), v.tmp);

            // Omega_D,phi
            routines::dot( v.p.mu[1]*v.p.tau[1], v.f.gradN(1), v.gradPsip, 0., v.tmp3[0]);

            // Divergence
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    /// ----------------------with electron density --------------------///
    {"oexbe", "ExB vorticity term with electron density", false,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDot( v.p.mu[1], result, v.f.density(0), 0., result);
        }
    },
    {"odiae", "Diamagnetic vorticity term with electron density", false,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( v.f.gradN(0), v.gradPsip, result);
            dg::blas1::scal( result, v.p.mu[1]*v.p.tau[1]);
        }
    },
    {"jsoexbe_tt", "ExB vorticity flux term with electron density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB Dot GradPsi
            routines::jacobian( v.f.bhatgB(), v.f.gradP(0), v.gradPsip, v.tmp[0]);

            // Omega_E
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDot( v.p.mu[1], result, v.f.density(0), 0., result);

            // Multiply everything
            dg::blas1::pointwiseDot( 1., result, v.tmp[0], 0., result);
        }
    },
    {"divoexbe_tt", "ExB vorticity flux term with electron density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB Dot GradPsi
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp); //u_E

            // Omega_E
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDot( v.p.mu[1], result, v.f.density(0), 0., v.tmp3[0]);

            // Divergence
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"jsoexbeUD_tt", "ExB vorticity flux term by diamagnetic velocity with electron density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // bxGradN/B Dot GradPsi
            routines::jacobian( v.p.tau[1], v.f.bhatgB(), v.f.gradN(0), v.gradPsip, 0., result);

            // m Omega_E,phi
            routines::dot( v.f.gradP(0), v.gradPsip, v.tmp[0]);
            dg::blas1::pointwiseDot( v.p.mu[1], v.tmp[0], v.f.binv(), v.f.binv(), 0., v.tmp[0]);

            // Multiply everything
            dg::blas1::pointwiseDot( 1., result, v.tmp[0], 0., result);
        }
    },
    {"divoexbeUD_tt", "ExB vorticity flux term by diamagnetic velocity with electron density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            routines::times( v.f.bhatgB(), v.f.gradN(0), v.tmp);

            // Omega_E
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( v.p.tau[1]*v.p.mu[1], result, v.f.binv(), v.f.binv(), 0., v.tmp3[0]);

            // Divergence
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"jsodiaeUE_tt", "Diamagnetic vorticity flux by ExB velocity term with electron density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB Dot GradPsi
            routines::jacobian( v.f.bhatgB(), v.f.gradP(0), v.gradPsip, result);

            // Omega_D,phi
            routines::dot( v.p.mu[1]*v.p.tau[1], v.f.gradN(0), v.gradPsip, 0., v.tmp[0]);

            // Multiply everything
            dg::blas1::pointwiseDot( 1., result, v.tmp[0], 0., result);
        }
    },
    {"divodiaeUE_tt", "Diamagnetic vorticity flux by ExB veloctiy term with electron density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB
            routines::times( v.f.bhatgB(), v.f.gradP(0), v.tmp);

            // Omega_D,phi
            routines::dot( v.p.mu[1]*v.p.tau[1], v.f.gradN(0), v.gradPsip, 0., v.tmp2[1]);

            // Divergence
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    /// ----------------------Remainders--------------------------------///
    {"jsoApar_tt", "A parallel vorticity flux term (Maxwell stress) (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            if( v.p.beta == 0) // avoid division by zero
            {
                dg::blas1::scal( result, 0.);
            }
            else
            {
                routines::jacobian( v.f.bhatgB(), v.f.gradA(), v.gradPsip, result);
                routines::dot( -1./v.p.beta, v.f.gradA(), v.gradPsip, 0., v.tmp[0]);
                dg::blas1::pointwiseDot( result, v.tmp[0], result);
            }
        }
    },
    {"divoApar_tt", "A parallel vorticity flux term (Maxwell stress) (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            if( v.p.beta == 0) // avoid division by zero
            {
                dg::blas1::scal( result, 0.);
            }
            else
            {
                routines::times( v.f.bhatgB(), v.f.gradA(), v.tmp);
                routines::dot( -1./v.p.beta, v.f.gradA(), v.gradPsip, 0., v.tmp3[0]);
                v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
            }
        }
    },
    {"jsodiaApar_tt", "A parallel diamagnetic vorticity flux term (magnetization stress) (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( v.gradPsip, v.f.gradU(1), v.tmp[0]);
            routines::dot( v.gradPsip, v.f.gradN(1), v.tmp[1]);
            dg::blas1::pointwiseDot( 1., v.tmp[0], v.f.density(1), 1., v.tmp[0], v.f.velocity(1), 0., result);

            routines::jacobian( v.f.bhatgB(), v.f.gradA(), v.gradPsip, result);
            dg::blas1::pointwiseDot( -1./2.*v.p.tau[1], result, v.tmp[0], 0., result);
        }
    },
    {"jsoexbApar_tt", "A parallel ExB vorticity flux term (magnetization stress) (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            routines::jacobian( v.f.bhatgB(), v.f.gradU(1), v.gradPsip, v.tmp[0]);
            routines::jacobian( v.f.bhatgB(), v.f.gradN(1), v.gradPsip, v.tmp[1]);
            dg::blas1::pointwiseDot( 1., v.tmp[0], v.f.density(1), 1., v.tmp[1], v.f.velocity(1), 0., result);
            routines::dot( v.f.gradA(), v.gradPsip, v.tmp[2]);
            dg::blas1::pointwiseDot( -1./2.*v.p.tau[1], result, v.tmp[2], 0., result);
        }
    },
    {"socurve_tt", "Vorticity source term electron curvature (Time average)", true,
        []( dg::x::DVec& result, Variables& v) {
            routines::dot( v.f.curv(), v.gradPsip, result);
            dg::blas1::pointwiseDot( -v.p.tau[0], v.f.density(0), result, 0., result);
        }
    },
    {"socurvi_tt", "Vorticity source term ion curvature (Time average)", true,
        []( dg::x::DVec& result, Variables& v) {
            routines::dot( v.f.curv(), v.gradPsip, result);
            dg::blas1::pointwiseDot( v.p.tau[1], v.f.density(1), result, 0., result);
        }
    },
    {"socurvkappae_tt", "Vorticity source term electron kappa curvature (Time average)", true,
        []( dg::x::DVec& result, Variables& v) {
            routines::dot( v.f.curvKappa(), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., v.f.density(0), v.f.velocity(0), v.f.velocity(0), 0., v.tmp[0]);
            dg::blas1::pointwiseDot( -v.p.mu[0], v.tmp[0], result, 0., result);
        }
    },
    {"socurvkappai_tt", "Vorticity source term ion kappa curvature (Time average)", true,
        []( dg::x::DVec& result, Variables& v) {
            routines::dot( v.f.curvKappa(), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., v.f.density(1), v.f.velocity(1), v.f.velocity(1), 0., v.tmp[0]);
            dg::blas1::pointwiseDot( v.p.mu[1], v.tmp[0], result, 0., result);
        }
    },
    {"sosne_tt", "ExB vorticity source term with electron source", true,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDot( v.p.mu[1], result, v.f.density_source(0), 0., result);
        }
    },
    {"sospi_tt", "Diamagnetic vorticity source term with electron source", true,
        []( dg::x::DVec& result, Variables& v){
            v.f.compute_gradSN( 0, v.tmp);
            routines::dot( v.p.mu[1]*v.p.tau[1], v.tmp, v.gradPsip, 0., result);
        }
    },
    {"loexbe_tt", "Vorticity dissipation term with electron Lambda", true,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( v.f.gradP(0), v.gradPsip, result);
            dg::blas1::pointwiseDot( 1., result, v.f.binv(), v.f.binv(), 0., result);

            v.f.compute_perp_diffusiveN( 1., v.f.density(0), v.tmp[0],
                    v.tmp[1], 0., v.tmp[2]);
            dg::blas1::axpby( v.p.nu_parallel_n, v.f.lapParN(0), 0., v.tmp[1]);
            dg::blas1::pointwiseDot( v.p.mu[1], v.tmp[1], result,
                v.p.mu[1], v.tmp[2], result, 0., result);
        }
    }
};

std::vector<Record> ParallelMomDiagnostics2d_list = { //36
    ///-----------------------Parallel momentum terms ------------------------//
    {"niui", "Product of ion gyrocentre density and velocity", false,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot(
                v.f.density(1), v.f.velocity(1), result);
        }
    },
    {"jsparexbi_tt", "Parallel momentum radial flux by ExB velocity with electron potential (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB Dot GradPsi
            routines::jacobian( v.f.bhatgB(), v.f.gradP(0), v.gradPsip, result);

            // parallel momentum mu_iN_iU_i
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1), v.f.velocity(1), result, 0., result);
        }
    },
    {"divparexbi_tt", "Divergence of parallel momentum radial flux by ExB velocity with electron potential (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp); //u_E
            // parallel momentum mu_iN_iU_i
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1), v.f.velocity(1), 0., v.tmp3[0]);
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"divparexbii_tt", "Divergence of parallel momentum radial flux by ExB velocity with ion potential (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB
            routines::times(v.f.bhatgB(), v.f.gradP(1), v.tmp); //u_E
            // parallel momentum mu_iN_iU_i
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1), v.f.velocity(1), 0., v.tmp3[0]);
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"divpardiai_tt", "Parallel momentum radial flux by Diamagnetic velocity with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( 1., v.f.curv(), v.f.gradN(1), 0., v.tmp[0]);
            routines::dot( 1., v.f.curv(), v.f.gradU(1), 0., v.tmp[1]);
            // Multiply everything
            dg::blas1::pointwiseDot( v.p.mu[1]*v.p.tau[1], v.tmp[0], v.f.velocity(1),
                v.p.mu[1]*v.p.tau[1], v.tmp[1], v.f.density(1), 0., result);
        }
    },
    {"divparkappai_tt", "Parallel momentum radial flux by curvature velocity (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1), v.f.velocity(1), 0., v.tmp[0]); // mu NU
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.velocity(1), v.f.velocity(1), v.tmp[0], 0., v.tmp[1]); //muNU mu U**2
            // mu NU(mu U^2 + 2tau)
            dg::blas1::axpbypgz( 2.*v.p.tau[1], v.tmp[0], +1., v.tmp[1], 0., v.tmp3[0]);
            v.f.centered_div( v.tmp3[0], v.f.curvKappa(), v.tmp2[0], result);
        }
    },
    {"divparmirrorAi_tt", "Divergence of parallel magnetic flutter force (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            //b_\perp
            v.f.compute_bperp( v.tmp);
            v.f.centered_div( v.f.density(1), v.tmp, v.tmp2[0], result);
            dg::blas1::scal( result, v.p.tau[1]);
        }
    },
    {"divparmirrorAe_tt", "Divergence of parallel magnetic flutter force (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            //b_\perp
            v.f.compute_bperp( v.tmp);
            v.f.centered_div( v.f.density(0), v.tmp, v.tmp2[0], result);
            dg::blas1::scal( result, -1.);
        }
    },
    {"divparApari_tt", "Parallel momentum radial flux by magnetic flutter (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.velocity(1), v.f.velocity(1), v.f.density(1),  0., v.tmp3[0]);
            //b_\perp^v
            v.f.compute_bperp( v.tmp);
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"divparApare_tt", "Divergence of parallel momentum radial flux by magnetic flutter (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::pointwiseDot( -v.p.mu[0], v.f.velocity(0), v.f.velocity(0), v.f.density(0),  0., v.tmp3[0]);
            //b_\perp^v
            v.f.compute_bperp( v.tmp);
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    /// --------------------- Parallel momentum source terms ---------------------//
    {"divjpari_tt", "Divergence of parallel ion momentum flux", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.divNUb(1), v.f.velocity(1),
                    0., result);
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1),
                    v.f.velocity(1), v.f.dsU(1), 1., result);
        }
    },
    {"divjpare_tt", "Divergence of parallel electron momentum flux", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot( -v.p.mu[0], v.f.divNUb(0), v.f.velocity(0),
                    0., result);
            dg::blas1::pointwiseDot( -v.p.mu[0], v.f.density(0),
                    v.f.velocity(0), v.f.dsU(0), 1., result);
        }
    },
    {"lparpar_tt", "Parallel momentum dissipation by parallel diffusion", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_parallel_diffusiveN( 1, v.tmp[0]);
            v.f.compute_parallel_diffusiveU( 1, v.tmp[1]);
            dg::blas1::pointwiseDot( 1., v.f.density(1), v.tmp[1], 1.,
                v.f.velocity(1), v.tmp[0], 0., result);
        }
    },
    {"lparperp_tt", "Parallel momentum dissipation by perp diffusion", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_perp_diffusiveN( 1., v.f.density(1), v.tmp3[0], v.tmp3[1],
                    0., v.tmp[0]);
            v.f.compute_perp_diffusiveU( 1., v.f.velocity(1), v.f.density(1),
                    v.tmp3[0], v.tmp3[1], v.tmp3[2], v.tmp2[0], 0., v.tmp[1]);
            dg::blas1::pointwiseDot( 1., v.tmp[0], v.f.velocity(1),
                    1., v.tmp[1], v.f.density(1), 0., result);
        }
    },
    {"lparparbphi_tt", "Parallel angular momentum dissipation by parallel diffusion", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_parallel_diffusiveN( 1, v.tmp[0]);
            v.f.compute_parallel_diffusiveU( 1, v.tmp[1]);
            dg::blas1::pointwiseDot( 1., v.f.density(1), v.tmp[1], 1.,
                v.f.velocity(1), v.tmp[0], 0., result);
            dg::blas1::pointwiseDot( result, v.f.bphi(), result);
        }
    },
    {"lparperpbphi_tt", "Parallel angular momentum dissipation by perp diffusion", true,
        []( dg::x::DVec& result, Variables& v ) {
            v.f.compute_perp_diffusiveN( 1., v.f.density(1), v.tmp3[0], v.tmp3[1],
                    0., v.tmp[0]);
            v.f.compute_perp_diffusiveU( 1., v.f.velocity(1), v.f.density(1),
                    v.tmp3[0], v.tmp3[1], v.tmp3[2], v.tmp2[0], 0., v.tmp[1]);
            dg::blas1::pointwiseDot( 1., v.tmp[0], v.f.velocity(1),
                    1., v.tmp[1], v.f.density(1), 0., result);
            dg::blas1::pointwiseDot( result, v.f.bphi(), result);
        }
    },
    //not so important
    {"sparKappaphii_tt", "Kappa Phi Source for parallel momentum", true,
        []( dg::x::DVec& result, Variables& v ) {
            routines::dot( v.f.curvKappa(), v.f.gradP(1), result);
            dg::blas1::pointwiseDot( -v.p.mu[1], v.f.density(1), v.f.velocity(1), result, 0., result);
        }
    },
    //not so important
    {"sparmirrorKappai_tt", "Generalized mirror force (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::pointwiseDot( v.p.mu[1]*v.p.tau[1], v.f.density(1), v.f.velocity(1),
                v.f.divCurvKappa(), 0., result);
        }
    },
    ///-----------------------Parallel anbular momentum terms-----------------//
    {"niuibphi", "Product of NiUi and covariant phi component of magnetic field unit vector", false,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot( 1.,
                v.f.density(1), v.f.velocity(1), v.f.bphi(), 0., result);
        }
    },
    {"jsparbphiexbi_tt", "Parallel angular momentum radial flux by ExB velocity with electron potential (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB Dot GradPsi
            routines::jacobian( v.f.bhatgB(), v.f.gradP(0), v.gradPsip, result);

            // parallel momentum mu_iN_iU_i
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1), v.f.velocity(1), result, 0., result);

            // Multiply bphi
            dg::blas1::pointwiseDot( 1., result, v.f.bphi(), 0., result);
        }
    },
    {"divparbphiexbi_tt", "Divergence of parallel angular momentum radial flux by ExB velocity with electron potential (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp); //u_E
            // parallel momentum mu_iN_iU_i bphi
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1), v.f.velocity(1), v.f.bphi(), 0., v.tmp3[0]);
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"divparbphiexbii_tt", "Divergence of parallel angular momentum radial flux by ExB velocity with ion potential (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            // ExB
            routines::times(v.f.bhatgB(), v.f.gradP(1), v.tmp); //u_E
            // parallel momentum mu_iN_iU_i bphi
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1), v.f.velocity(1), v.f.bphi(), 0., v.tmp3[0]);
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"divparbphidiai_tt", "Parallel angular momentum radial flux by Diamagnetic velocity with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            routines::dot( 1., v.f.curv(), v.f.gradN(1), 0., v.tmp[0]);
            routines::dot( 1., v.f.curv(), v.f.gradU(1), 0., v.tmp[1]);
            dg::blas1::pointwiseDot( v.tmp[0], v.f.bphi(), v.tmp[0]);
            dg::blas1::pointwiseDot( v.tmp[1], v.f.bphi(), v.tmp[1]);
            // Multiply everything
            dg::blas1::pointwiseDot( v.p.mu[1]*v.p.tau[1], v.tmp[0], v.f.velocity(1),
                v.p.mu[1]*v.p.tau[1], v.tmp[1], v.f.density(1), 0., result);
        }
    },
    {"divparbphikappai_tt", "Parallel angular momentum radial flux by curvature velocity (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.density(1), v.f.velocity(1), 0., v.tmp[0]);
            dg::blas1::pointwiseDot( v.p.mu[1], v.f.velocity(1), v.f.velocity(1), v.tmp[0],  0., v.tmp[1]);
            // mu NU(mu U^2 + 2tau)bphi
            dg::blas1::pointwiseDot( 2.*v.p.tau[1], v.tmp[0], v.f.bphi(), +1., v.tmp[1], v.f.bphi(), 0., v.tmp3[0]);
            v.f.centered_div( v.tmp3[0], v.f.curvKappa(), v.tmp2[0], result);
        }
    },
    {"divparbphiApar_tt", "Parallel angular momentum radial flux by magnetic flutter (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::pointwiseDot( -v.p.mu[0], v.f.velocity(0), v.f.velocity(0), v.f.density(0),  0., result);
            dg::blas1::pointwiseDot( +v.p.mu[1], v.f.velocity(1), v.f.velocity(1), v.f.density(1),  1., result);
            dg::blas1::axpbypgz( -v.p.tau[0], v.f.density(0),
                                 +v.p.tau[1], v.f.density(1), 1., result);
            dg::blas1::pointwiseDot( v.f.bphi(), result, v.tmp3[0]);
            v.f.compute_bperp( v.tmp);
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    /// --------------------electron force balance usually well-fulfilled ----//
    {"sparphie_tt", "Electric force in electron momentum density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::pointwiseDot( 1., v.f.dsP(0), v.f.density(0), 0., result);
        }
    },
    {"sparmirrore_tt", "Parallel electron pressure (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::axpby( v.p.tau[0], v.f.dsN(0), 0., result);
        }
    },
    {"sparmirrorAe_tt", "Apar Mirror force term with electron density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            v.f.compute_bperp( v.tmp);
            routines::dot( v.p.tau[0], v.tmp, v.f.gradN(0), 0., result);
        }
    },

    {"sparphiAe_tt", "Apar Electric force in electron momentum density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            v.f.compute_bperp( v.tmp);
            routines::dot( v.tmp, v.f.gradP(0), result);
            dg::blas1::pointwiseDot( v.f.density(0), result, result);
        }
    },
    {"spardotAe_tt", "Apar Electric force in electron momentum density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            v.f.compute_dot_aparallel( result);
            dg::blas1::pointwiseDot( v.f.density(0), result, result);
        }
    },
    {"neue", "Product of electron density and velocity", false,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot(
                v.f.density(0), v.f.velocity(0), result);
        }
    },
    /// -----------Ion force balance ----------------------///
    {"sparphii_tt", "Electric force term in ion momentum density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::pointwiseDot( -1., v.f.dsP(1), v.f.density(1), 0., result);
        }
    },
    {"sparmirrori_tt", "Parallel ion pressure (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            dg::blas1::axpby( -v.p.tau[1], v.f.dsN(1), 0., result);
        }
    },
    {"sparmirrorAi_tt", "Apar Mirror force term with ion density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            v.f.compute_bperp( v.tmp);
            routines::dot( -v.p.tau[1], v.tmp, v.f.gradN(1), 0., result);
        }
    },
    {"sparphiAi_tt", "Apar Electric force in ion momentum density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            v.f.compute_bperp( v.tmp);
            routines::dot( v.tmp, v.f.gradP(1), result);
            dg::blas1::pointwiseDot( -1., v.f.density(1), result, 0., result);
        }
    },
    {"spardotAi_tt", "Apar Electric force in ion momentum density (Time average)", true,
        []( dg::x::DVec& result, Variables& v){
            v.f.compute_dot_aparallel( result);
            dg::blas1::pointwiseDot( -1., v.f.density(1), result, 0., result);
        }
    },
    {"friction_tt", "Friction force in momentum density (Time average)", true,
        []( dg::x::DVec& result, Variables& v ) {
            dg::blas1::pointwiseDot( 1., v.f.velocity(1), v.f.density(1), -1.,
                    v.f.velocity(0), v.f.density(0), 0., result);
            dg::blas1::pointwiseDot( v.p.eta, result, v.f.density(0), 0, result);
        }
    },
};

std::vector<Record> RSDiagnostics2d_list = { //2
    /// --------------------- Zonal flow energy terms------------------------//
    {"nei0", "inertial factor", false,
        []( dg::x::DVec& result, Variables& v ) {
            routines::dot( v.gradPsip, v.gradPsip, result);
            dg::blas1::transform( result, result, feltor::Positive());
            dg::blas1::pointwiseDot( v.mag.R0()*v.mag.R0(),
                result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDivide( v.f.density(0), result, result);
        }
    },
    {"snei0_tt", "inertial factor source", true,
        []( dg::x::DVec& result, Variables& v ) {
            routines::dot( v.gradPsip, v.gradPsip, result);
            dg::blas1::transform( result, result, feltor::Positive());
            dg::blas1::pointwiseDot( v.mag.R0()*v.mag.R0(),
                result, v.f.binv(), v.f.binv(), 0., result);
            dg::blas1::pointwiseDivide( v.f.density_source(0), result, result);
        }
    }
};

std::vector<Record> COCEDiagnostics2d_list = { // 16
    /// ----------------- COCE EQUATION ----------------//
    /// ---------- Polarization charge densities -----------///
    {"v_Omega_E", "Electron polarisation term", false,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_pol( 1., v.f.density(0), v.tmp[0], 0., result);
        }
    },
    {"v_Omega_E_gf", "Ion polarisation term", false,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_pol( 1., v.f.density(1), v.tmp[0], 0., result);
        }
    },
    /// ------------ Polarization advections ------------------//
    //The fsa of the main and rest terms is almost the same as the one of divoexbi
    {"v_adv_E_main_tt", "Main electric advective term (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_pol( 1., v.f.density(0), v.tmp[0], 0., v.tmp3[0]);
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp); //u_E
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"v_adv_E_main_gf_tt", "Main electric advective term (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_pol( 1., v.f.density(1), v.tmp[0], 0., v.tmp3[0]);
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp); //u_E
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
        }
    },
    {"v_adv_E_rest_tt", "Electric advective term (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            // NOT implemented for true curvature mode
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp2); //u_E
            dg::blas1::pointwiseDot(1., v.f.binv(), v.f.binv(), v.f.density(0), 0., v.tmp[0]);
            routines::scal(v.tmp[0], v.f.gradP(0), v.tmp3); //ne Grad_phi/B^2
            v.f.centered_v_dot_nabla(v.tmp3, v.tmp2[0], v.tmp[2], v.tmp[0]); //t3*nabla(u_E^R)
            v.f.centered_v_dot_nabla(v.tmp3, v.tmp2[1], v.tmp[2], v.tmp[1]); //t3*nabla(u_E^Z)
            v.f.centered_div( v.p.mu[1], v.tmp, v.tmp2[0], result);
        }
    },
    {"v_adv_E_rest_gf_tt", "Electric advective term GF (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            // NOT implemented for true curvature mode
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp2); //u_E
            dg::blas1::pointwiseDot(1., v.f.binv(), v.f.binv(), v.f.density(1), 0., v.tmp[0]);
            routines::scal(v.tmp[0], v.f.gradP(0), v.tmp3); //ne Grad_phi/B^2
            v.f.centered_v_dot_nabla(v.tmp3, v.tmp2[0], v.tmp[2], v.tmp[0]); //t3*nabla(u_E^R)
            v.f.centered_v_dot_nabla(v.tmp3, v.tmp2[1], v.tmp[2], v.tmp[1]); //t3*nabla(u_E^Z)
            v.f.centered_div( v.p.mu[1], v.tmp, v.tmp2[0], result);
        }
    },
    //The fsa of the main and rest terms is almost the same as the one of divodiaiUE
    {"v_adv_D_main_tt", "Main diamagnetic term (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_lapMperpN(-1.0, v.f.density(0), v.tmp[0], 0., v.tmp3[0]);
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp); //u_E
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
            dg::blas1::scal( result, v.p.mu[1]*v.p.tau[1]);
        }
    },
    {"v_adv_D_main_gf_tt", "Main diamagnetic term (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_lapMperpN(-1.0, v.f.density(1), v.tmp[0], 0., v.tmp3[0]);
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp); //u_E
            v.f.centered_div( v.tmp3[0], v.tmp, v.tmp2[0], result);
            dg::blas1::scal( result, v.p.mu[1]*v.p.tau[1]);
        }
    },
    {"v_adv_D_rest_tt", "Diamagnetic advective term (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            // NOT implemented for true curvature mode
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp2); //u_E
            v.f.centered_v_dot_nabla(v.f.gradN(0), v.tmp2[0], v.tmp[2], v.tmp[0]); //t3*nabla(u_E^R)
            v.f.centered_v_dot_nabla(v.f.gradN(0), v.tmp2[1], v.tmp[2], v.tmp[1]); //t3*nabla(u_E^Z)
            v.f.centered_div( v.p.tau[1]*v.p.mu[1], v.tmp, v.tmp2[0], result);
        }
    },

    {"v_adv_D_rest_gf_tt", "Diamagnetic advective term GF (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            // NOT implemented for true curvature mode
            routines::times(v.f.bhatgB(), v.f.gradP(0), v.tmp2); //u_E
            v.f.centered_v_dot_nabla(v.f.gradN(1), v.tmp2[0], v.tmp[2], v.tmp[0]); //t3*nabla(u_E^R)
            v.f.centered_v_dot_nabla(v.f.gradN(1), v.tmp2[1], v.tmp[2], v.tmp[1]); //t3*nabla(u_E^Z)
            v.f.centered_div( v.p.tau[1]*v.p.mu[1], v.tmp, v.tmp2[0], result);
        }
    },
    ///---------------- J_b_perp components --------//
    {"v_J_mag_tt", "Magnetization current term GF (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            // take care to use correct derivatives...
            v.f.compute_lapMperpN( -1., v.f.density(1), v.tmp[0], 0., v.tmp2[0]);
            dg::blas1::pointwiseDot( v.f.velocity(1), v.tmp2[0], v.tmp2[0]);
            v.f.compute_lapMperpU( 1, v.tmp2[1]);
            dg::blas1::pointwiseDot( v.f.density(1), v.tmp2[1], v.tmp2[1]);
            dg::tensor::scalar_product3d( 1., 1.,
                    v.f.gradN(1)[0], v.f.gradN(1)[1], v.f.gradN(1)[2],
                    v.f.projection(), 1.,
                    v.f.gradU(1)[0], v.f.gradU(1)[1], v.f.gradU(1)[2],
                    0., v.tmp2[2]);
            v.f.compute_bperp(v.tmp);
            dg::blas1::axpbypgz( 1., v.tmp2[0], 1., v.tmp2[1], 1., v.tmp2[2]);
            v.f.centered_div( v.tmp2[2], v.tmp, v.tmp2[0], result);
            dg::blas1::scal ( result, v.p.tau[1]/2.);
        }
    },
    {"v_J_bperp_tt", "Div J_par times b_perp term (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_bperp(v.tmp);
            dg::blas1::pointwiseDot(1., v.f.density(1), v.f.velocity(1), -1.,
                    v.f.density(0), v.f.velocity(0), 0, v.tmp2[0]);
            v.f.centered_div( v.tmp2[0], v.tmp, v.tmp3[0], result);
        }
    },
    ///-------------- Sources term----------///
    {"v_S_E_tt", "Electric source vorticity (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_source_pol( 1., v.f.density_source(0), v.tmp[0], 0., result);
        }
    },
    {"v_S_D_tt", "Diamagnetic source vorticity (time integrated)", true,
        []( dg::x::DVec& result, Variables& v) {
            v.f.compute_lapMperpN(-v.p.tau[1]*v.p.mu[1], v.f.density_source(0), v.tmp[0], 0., result);
        }
    },
    /// ------------- Radial Force Balance terms ---------///
    {"RFB_E_r_GradPsip_tt", "Radial electric field in RFB (time integrated)", true,
         []( dg::x::DVec& result, Variables& v){
             routines::dot(-1., v.f.gradP(0), v.gradPsip, 0., result);
         }
     },
    {"RFB_GradPi_GradPsip_tt", "Radial pressure gradient component of RFB (time integrated)", true,
         []( dg::x::DVec& result, Variables& v){
             routines::dot(v.f.gradN(0), v.gradPsip, result);
             dg::blas1::pointwiseDivide(v.p.tau[1], result, v.f.density(0), 0.,
                     result);
         }
     },

};

// probes list
// Idea: I think we can unfold the species dependent lists into one long list
std::vector<dg::file::Record<void(dg::x::DVec&,Variables&)>> probe_list = {
     {true, "n", "probe measurement of density",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.density(0), result);
         }
     },
     {true, "u", "probe measurement of parallel velocity",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.velocity(0), result);
         }
     },
     {true, "phi", "probe measurement of electric potential",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.potential(0), result);
         }
     },
     {false, "apar", "probe measurement of parallel magnetic potential",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.apar(), result);
         }
     },
     {"n_R", "probe measurement of d/dR electron density",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradN(0)[0], result);
         }
     },
     {"u_R", "probe measurement of d/dR parallel electron velocity",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradU(0)[0], result);
         }
     },
     {"phi_R", "probe measurement of d/dR electric potential",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradP(0)[0], result);
         }
     },
     {"aparR", "probe measurement of d/dR parallel magnetic potential",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradA()[0], result);
         }
     },
     {"neZ", "probe measurement of d/dZ electron density",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradN(0)[1], result);
         }
     },
     {"niZ", "probe measurement of d/dZ ion density",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradN(1)[1], result);
         }
     },
     {"ueZ", "probe measurement of d/dZ parallel electron velocity",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradU(0)[1], result);
         }
     },
     {"uiZ", "probe measurement of d/dZ parallel ion velocity",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradU(1)[1], result);
         }
     },
     {"phiZ", "probe measurement of d/dZ electric potential",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradP(0)[1], result);
         }
     },
     {"aparZ", "probe measurement of d/dZ parallel magnetic potential",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.gradA()[1], result);
         }
     },
     {"nePar", "probe measurement of d/dPar electron density",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.dsN(0), result);
         }
     },
     {"niPar", "probe measurement of d/dPar ion density",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.dsN(1), result);
         }
     },
     {"uePar", "probe measurement of d/dPar parallel electron velocity",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.dsU(0), result);
         }
     },
     {"uiPar", "probe measurement of d/dPar parallel ion velocity",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.dsU(1), result);
         }
     },
     {"phiPar", "probe measurement of d/dPar electric potential",
         []( dg::x::DVec& result, Variables& v ) {
              dg::blas1::copy(v.f.dsP(0), result);
         }
     }
 };
*/

///%%%%%%%%%%%%%%%%%%%%%%%%%%END DIAGNOSTICS LIST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///%%%%%%%%%%%%%%%%%%%%%%%%%%END DIAGNOSTICS LIST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///%%%%%%%%%%%%%%%%%%%%%%%%%%END DIAGNOSTICS LIST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template<class NcFile>
void write_global_attributes(NcFile& file, int argc, char* argv[], std::string inputfile)
{
    // called only by master thread
    std::map<std::string, dg::file::nc_att_t> att;
    att["title"] = "Output file of feltor/src/thermal/thermal.cpp";
    att["Conventions"] = "CF-1.8";
    att["history"] = dg::file::timestamp( argc, argv);
    att["comment"] = "Find more info in feltor/src/thermal/thermal.tex";
    att["source"] = "FELTOR";
    att["references"] = "https://github.com/feltor-dev/feltor";
    att["inputfile"] = inputfile;
    file.put_atts( att);

    file.put_atts( dg::file::version_flags);
}

template<class NcFile, class HostList, class HostList2>
void write_static_list( NcFile& file, const HostList& records, const HostList2& species_records,
    Variables& var,
    const dg::x::CylindricalGrid3d& grid, const dg::x::CylindricalGrid3d& g3d_out,
    dg::geo::CylindricalFunctor transition, const std::vector<std::string>& species_names )
{
    // the unique thing here is that we evaluate 3d but only write 2d
    std::unique_ptr<dg::x::aGeometry2d> g2d_out_ptr( g3d_out.perp_grid());

    dg::x::HVec resultH = dg::evaluate( dg::zero, grid);
    dg::x::HVec transferH( dg::evaluate(dg::zero, g3d_out));
    dg::MultiMatrix<dg::x::HMatrix,dg::x::HVec> projectH =
            dg::create::fast_projection( grid, 1, var.p.cx, var.p.cy);
    for ( auto& record : records)
    {
        record.function( resultH, var, grid);
        dg::blas2::symv( projectH, resultH, transferH);
        file.defput_var( record.name, {"y", "x"}, record.atts, {*g2d_out_ptr},
                transferH);
    }
    resultH = dg::pullback( transition, grid);
    dg::blas2::symv( projectH, resultH, transferH);
    file.defput_var( "MagneticTransition", {"y", "x"}, {{"long_name",
        "The region where the magnetic field is modified"}},
        {*g2d_out_ptr}, transferH);
    for( unsigned s=0; s<species_names.size(); s++)
    for( auto& record: species_records)
    {
        record.function( resultH, var, grid, s);
        dg::blas2::symv( projectH, resultH, transferH);
        file.defput_var( species_names[s] + "_" + record.name, {"y", "x"}, record.atts,
                {*g2d_out_ptr}, transferH);
    }
}

void append_equations( std::vector<thermal::Record>& list, const std::vector<thermal::Record>& b)
{
    list.insert( list.begin(), b.begin(), b.end());
}
std::vector<thermal::Record> generate_equation_list( const dg::file::WrappedJsonValue& js)
{
    std::vector<thermal::Record> list;
    bool equation_list_exists = js["output"].isMember("equations");
    if( equation_list_exists)
    {
        for( unsigned i=0; i<js["output"]["equations"].size(); i++)
        {
            std::string eqn = js["output"]["equations"][i].asString();
            if( eqn == "Basic")
                append_equations( list, thermal::basicDiagnostics2d_list);
            else if( eqn == "Mass-conserv")
                ;// append_equations( list, thermal::MassConsDiagnostics2d_list);
            else if( eqn == "Energy-theorem")
                ;// append_equations( list, thermal::EnergyDiagnostics2d_list);
            else if( eqn == "Toroidal-momentum")
                ;// append_equations( list, thermal::ToroidalExBDiagnostics2d_list);
            else if( eqn == "Parallel-momentum")
                ;// append_equations( list, thermal::ParallelMomDiagnostics2d_list);
            else if( eqn == "Zonal-Flow-Energy")
                ;// append_equations( list, thermal::RSDiagnostics2d_list);
            else if( eqn == "COCE")
                ;// append_equations( list, thermal::COCEDiagnostics2d_list);
            else
                throw std::runtime_error( "output: equations: "+eqn+" not recognized!\n");
        }
    }
    else // default diagnostics
    {
        append_equations(list, thermal::basicDiagnostics2d_list);
        //append_equations(list, thermal::MassConsDiagnostics2d_list);
        //append_equations(list, thermal::EnergyDiagnostics2d_list);
        //append_equations(list, thermal::ToroidalExBDiagnostics2d_list);
        //append_equations(list, thermal::ParallelMomDiagnostics2d_list);
        //append_equations(list, thermal::RSDiagnostics2d_list);
    }
    return list;
}

// TODO Same as feltor::WriteIntegrateDiagnostics2dList except capabilities to loop multiple species
// Maybe we should rather unfold records list and then use feltor::Write ?
template< class NcFile>
struct WriteIntegrateDiagnostics2dList
{
    WriteIntegrateDiagnostics2dList( NcFile& file,
        const dg::x::CylindricalGrid3d& grid,
        const dg::x::CylindricalGrid3d& g3d_out,
        const std::vector<thermal::Record>& equation_list,
        const std::vector<std::string>& species_names) :
            m_file(&file), m_slab(grid), m_grid(grid), m_g3d_out(g3d_out),
            m_equation_list(equation_list), m_species_names(species_names)
    {
#ifdef WITH_MPI
        int rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif //WITH_MPI
        std::unique_ptr<dg::x::aGeometry2d> g2d_out_ptr( g3d_out.perp_grid());
        for( auto& record : m_equation_list)
        for( unsigned s=0; s<species_names.size(); s++)
        {
            std::string name = record.name + "_ta2d";
            std::string long_name = record.long_name + " (Toroidal average)";
            if( record.species_dependent)
                m_file->template def_var_as<double>( species_names[s] + "_" + name,
                    {"time", "y", "x"}, {{"long_name", long_name}});
            else if( s == 0)
                m_file->template def_var_as<double>( name,
                    {"time", "y", "x"}, {{"long_name", long_name}});
            name = record.name + "_2d";
            long_name = record.long_name+ " (Evaluated on phi = 0 plane)";
            if( record.species_dependent)
                m_file->template def_var_as<double>( species_names[s] + "_" + name,
                    {"time", "y", "x"}, {{"long_name", long_name}});
            else if( s == 0)
                m_file->template def_var_as<double>( name,
                    {"time", "y", "x"}, {{"long_name", long_name}});
        }
        m_slab = {*g2d_out_ptr};
        #ifdef WITH_MPI // only root group needs to track
        if( dg::file::detail::mpi_comm_global2local_rank( g2d_out_ptr->communicator()) == MPI_UNDEFINED)
            m_track  = false;
        #endif
        m_resultD = dg::evaluate( dg::zero, grid);
        m_transferD = dg::evaluate(dg::zero, g3d_out);
        m_transferH = dg::evaluate(dg::zero, g3d_out);
        m_projectD = dg::create::fast_projection( grid, 1, grid.Nx()/g3d_out.Nx(), grid.Ny()/g3d_out.Ny());
        m_transferH2d = dg::evaluate( dg::zero, *g2d_out_ptr);
        m_toroidal_average = { g3d_out, dg::coo3d::z};
    }
    // same as buffer and flush
    void write( double time, Variables& var)
    {
        buffer( time, var);
        flush(var);
    }

    void buffer( double time, Variables& var)
    {
        // evaluates function and updates time integrals for all integrals
        auto transferD2d_view = dg::split( m_transferD, m_g3d_out);
        for( auto& record : m_equation_list)
        for( unsigned s=0; s<m_species_names.size(); s++)
        {
            if( s!=0 and not record.species_dependent)
                break;
            if( record.integral)
            {
                record.function( m_resultD, var, s);
                dg::blas2::symv( m_projectD, m_resultD, m_transferD);
                //toroidal average and add to time integral
                std::string name = record.name+"_ta2d";
                if( record.species_dependent)
                    name = m_species_names[s] + "_"  + name;
                dg::assign( m_transferD, m_transferH);
                m_toroidal_average( m_transferH, m_transferH2d, false);
                if(m_track && m_first_buffer) m_time_integrals[name].init(
                    time, m_transferH2d);
                if(m_track && !m_first_buffer) m_time_integrals.at(name).add( time,
                    m_transferH2d);

                // 2d data of plane varphi = 0
                name = record.name + "_2d";
                if( record.species_dependent)
                    name = m_species_names[s] + "_"  + name;
                dg::split( m_transferD, transferD2d_view, m_g3d_out);
                dg::assign( transferD2d_view[0], m_transferH2d);
                if(m_track && m_first_buffer) m_time_integrals[name].init(
                    time, m_transferH2d);
                if(m_track && !m_first_buffer) m_time_integrals.at(name).add( time,
                    m_transferH2d);
            }
        }
        if( m_first_buffer)
            m_first_buffer = false;
    }
    void flush( Variables& var )
    {
        // write time integrals for
        auto transferD2d_view = dg::split( m_transferD, m_g3d_out);
        for( auto& record : m_equation_list)
        for( unsigned s=0; s<m_species_names.size(); s++)
        {
            if( s!=0 and not record.species_dependent)
                break;
            std::array<std::string, 2> names = { record.name+"_ta2d", record.name+"_2d"};
            if( record.species_dependent)
            {
                names[0] = m_species_names[s] + "_"  + names[0];
                names[1] = m_species_names[s] + "_"  + names[1];
            }
            if(record.integral) // we already computed the output...
            {
                for( std::string name : names)
                {
                    if(m_track) m_transferH2d = m_time_integrals.at(name).get_integral();
                    m_file->put_var( name, {m_start, m_slab}, m_transferH2d);
                    if(m_track) m_time_integrals.at(name).flush();
                }
            }
            else // compute from scratch
            {
                record.function( m_resultD, var, s);
                dg::blas2::symv( m_projectD, m_resultD, m_transferD);

                dg::assign( m_transferD, m_transferH);
                m_toroidal_average( m_transferH, m_transferH2d, false);
                m_file->put_var( names[0], {m_start, m_slab}, m_transferH2d);

                // 2d data of plane varphi = 0
                dg::split( m_transferD, transferD2d_view, m_g3d_out);
                dg::assign( transferD2d_view[0], m_transferH2d);
                m_file->put_var( names[1], {m_start, m_slab}, m_transferH2d);
            }
        }
        m_start++;
    }
    private:
    NcFile * m_file;
    typename NcFile::Hyperslab m_slab;
    bool m_first_buffer = true;
    bool m_track = true;
    size_t m_start = 0;
    dg::x::DVec m_resultD;
    dg::x::DVec m_transferD;
    dg::x::HVec m_transferH;
    dg::MultiMatrix<dg::x::DMatrix,dg::x::DVec> m_projectD;
    dg::x::HVec m_transferH2d;
    dg::Average<dg::x::IHMatrix, dg::x::HVec> m_toroidal_average;
    std::map<std::string, dg::Simpsons<dg::x::HVec>> m_time_integrals;
    const dg::x::CylindricalGrid3d m_grid, m_g3d_out;
    std::vector<thermal::Record> m_equation_list;
    std::vector<std::string> m_species_names;

};


}//namespace thermal
