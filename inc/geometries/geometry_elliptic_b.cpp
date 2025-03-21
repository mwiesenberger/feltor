#include <iostream>
#include <memory>
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "dg/algorithm.h"
#include "dg/file/file.h"

#include "flux.h"
#include "solovev.h"
#include "guenter.h"
#include "simple_orthogonal.h"
#include "curvilinear.h"
#ifdef WITH_MPI
#include "mpi_curvilinear.h"
#endif
#include "testfunctors.h"

int main(int argc, char**argv)
{
#ifdef WITH_MPI
    MPI_Init( &argc, &argv);
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
    unsigned n, Nx, Ny, Nz;
    double psi_0, psi_1;
#ifdef WITH_MPI
    MPI_Comm comm;
    dg::mpi_init3d( dg::DIR, dg::PER, dg::PER, n, Nx, Ny, Nz, comm);
    if(rank==0)std::cout << "Type psi_0 (-20) and psi_1 (-4)\n";
    if(rank==0)std::cin >> psi_0>> psi_1;
    MPI_Bcast( &psi_0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( &psi_1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    std::cout << "Type n (3), Nx (8), Ny (80), Nz (1)\n";
    std::cin >> n>> Nx>>Ny>>Nz;
    std::cout << "Type psi_0 (-20) and psi_1 (-4)\n";
    std::cin >> psi_0>> psi_1;
#endif
    auto js = dg::file::file2Json( argc == 1 ? "geometry_params_Xpoint.json" : argv[1]);
    //write parameters from file into variables
    dg::geo::solovev::Parameters gp(js);
    dg::geo::TokamakMagneticField mag = dg::geo::createSolovevField(gp);
    DG_RANK0 gp.display( std::cout);
    dg::Timer t;
    DG_RANK0 std::cout << "Psi min "<<mag.psip()(gp.R_0, 0)<<"\n";
    DG_RANK0 std::cout << "Constructing grid ... \n";
    t.tic();
    //dg::geo::SimpleOrthogonal generator( mag.get_psip(), psi_0, psi_1, gp.R_0, 0., 1);
    dg::geo::FluxGenerator generator( mag.get_psip(), mag.get_ipol(), psi_0, psi_1, gp.R_0, 0., 1);
    dg::geo::x::CurvilinearProductGrid3d g3d( generator, n, Nx, Ny,Nz, dg::DIR, dg::PER, dg::PER
#ifdef WITH_MPI
    , comm
#endif
    );
    std::unique_ptr<dg::x::aGeometry2d> g2d( g3d.perp_grid() );

    dg::Elliptic<dg::x::aGeometry2d, dg::x::DMatrix, dg::x::DVec> pol( *g2d,
        dg::forward);
    t.toc();
    DG_RANK0 std::cout << "Construction took "<<t.diff()<<"s\n";
    ///////////////////////////////////////////////////////////////////////////
    dg::file::NcFile file( "testE.nc", dg::file::nc_clobber);
    file.defput_dim( "x", {{"axis", "X"}}, g2d->abscissas(0));
    file.defput_dim( "y", {{"axis", "Y"}}, g2d->abscissas(1));
    file.defput_var( "xc", {"y", "x"}, {}, {*g2d}, g2d->map()[0]);
    file.defput_var( "yc", {"y", "x"}, {}, {*g2d}, g2d->map()[1]);

    ///////////////////////////////////////////////////////////////////////////
    dg::x::DVec x = dg::evaluate( dg::zero, *g2d);
    const dg::x::DVec b =    dg::pullback( dg::geo::EllipticDirPerM(mag, psi_0, psi_1, 4), *g2d);
    const dg::x::DVec chi =  dg::pullback( dg::geo::Bmodule(mag), *g2d);
    const dg::x::DVec solution = dg::pullback( dg::geo::FuncDirPer(mag, psi_0, psi_1, 4), *g2d);
    const dg::x::DVec vol3d = dg::create::volume( *g2d);
    pol.set_chi( chi);
    //compute error
    dg::x::DVec error( solution);
    const double eps = 1e-10;
    dg::PCG<dg::x::DVec > pcg( x, n*n*Nx*Ny*Nz);
    DG_RANK0 std::cout << "eps \t # iterations \t error \t hx_max\t hy_max \t time/iteration \n";
    DG_RANK0 std::cout << eps<<"\t";
    t.tic();
    unsigned number = pcg.solve(pol, x,b, pol.precond(), pol.weights(), eps);
    DG_RANK0 std::cout <<number<<"\t";
    t.toc();
    dg::blas1::axpby( 1.,x,-1., solution, error);
    double err = dg::blas2::dot( vol3d, error);
    const double norm = dg::blas2::dot( vol3d, solution);
    DG_RANK0 std::cout << sqrt( err/norm) << "\t";

    dg::SparseTensor<dg::x::DVec> metric = g2d->metric();
    dg::x::DVec gyy = metric.value(1,1), gxx=metric.value(0,0), volume = dg::tensor::volume(metric);
    dg::blas1::transform( gxx, gxx, dg::SQRT<double>());
    dg::blas1::transform( gyy, gyy, dg::SQRT<double>());
    dg::blas1::pointwiseDot( gxx, volume, gxx);
    dg::blas1::pointwiseDot( gyy, volume, gyy);
    dg::blas1::scal( gxx, g2d->hx());
    dg::blas1::scal( gyy, g2d->hy());
    double maxx = dg::blas1::reduce( gxx, -1e308, thrust::maximum<double>());
    DG_RANK0 std::cout << maxx << "\t";
    double maxy = dg::blas1::reduce( gyy, -1e308, thrust::maximum<double>());
    DG_RANK0 std::cout << maxy << "\t";
    DG_RANK0 std::cout<<t.diff()/(double)number<<"s"<<std::endl;
    ///////////////////////////////////////////////////////////////////////
    DG_RANK0 std::cout << "TESTING VARIATION\n";
    pol.variation( x, x);
    const dg::x::DVec variation = dg::pullback( dg::geo::VariationDirPer( mag, psi_0, psi_1), *g2d);
    dg::blas1::axpby( 1., variation, -1., x);
    double result = dg::blas2::dot( x, vol3d, x);
    DG_RANK0 std::cout << "               distance to solution "<<sqrt( result)<<std::endl; //don't forget sqrt when comuting errors

    file.defput_var( "error", {"y", "x"}, {}, {*g2d}, error);
    file.defput_var( "num_solution", {"y", "x"}, {}, {*g2d}, x);
    file.defput_var( "ana_solution", {"y", "x"}, {}, {*g2d}, solution);
    file.close();

#ifdef WITH_MPI
    MPI_Finalize();
#endif

    return 0;
}
