#pragma once

#include "exblas/exdot_omp.h"
#include "config.h"
#include "vector_categories.h"
#ifdef MPI_VERSION
#include "exblas/mpi_accumulate.h"
#endif //MPI_VERSION

namespace dg
{
template<class value_type>
void transpose_dispatch( OmpTag, unsigned nx, unsigned ny, const value_type* RESTRICT in, value_type* RESTRICT out)
{
#pragma omp parallel for
    for( unsigned i=0; i<ny; i++)
        for( unsigned j=0; j<nx; j++)
            out[j*ny+i] = in[i*nx+j];
}
template<class value_type>
void extend_line( OmpTag, unsigned nx, unsigned ny, const value_type* RESTRICT in, value_type* RESTRICT out)
{
#pragma omp parallel for
    for( unsigned i=0; i<ny; i++)
        for( unsigned j=0; j<nx; j++)
            out[i*nx+j] = in[j];
}
template<class value_type>
void extend_column( OmpTag, unsigned nx, unsigned ny, const value_type* RESTRICT in, value_type* RESTRICT out)
{
#pragma omp parallel for
    for( unsigned i=0; i<ny; i++)
        for( unsigned j=0; j<nx; j++)
            out[i*nx+j] = in[i];
}

template<class value_type>
void average( OmpTag, unsigned nx, unsigned ny, const value_type* in0, const value_type* in1, value_type* out)
{
    static_assert( std::is_same<value_type, double>::value, "Value type must be double!");
    static thrust::host_vector<int64_t> h_accumulator;
    h_accumulator.resize( ny*exblas::BIN_COUNT);
    for( unsigned i=0; i<ny; i++)
        exblas::exdot_omp(nx, &in0[i*nx], &in1[i*nx], &h_accumulator[i*exblas::BIN_COUNT]);
    for( unsigned i=0; i<ny; i++)
        out[i] = exblas::cpu::Round( &h_accumulator[i*exblas::BIN_COUNT]);
}

#ifdef MPI_VERSION
//local data plus communication
template<class value_type>
void average_mpi( OmpTag, unsigned nx, unsigned ny, const value_type* in0, const value_type* in1, value_type* out, MPI_Comm comm, MPI_Comm comm_mod, MPI_Comm comm_mod_reduce )
{
    static_assert( std::is_same<value_type, double>::value, "Value type must be double!");
    static thrust::host_vector<int64_t> h_accumulator;
    static thrust::host_vector<int64_t> h_accumulator2;
    h_accumulator2.resize( ny*exblas::BIN_COUNT);
    for( unsigned i=0; i<ny; i++)
        exblas::exdot_omp(nx, &in0[i*nx], &in1[i*nx], &h_accumulator2[i*exblas::BIN_COUNT]);
    h_accumulator.resize( h_accumulator2.size());
    exblas::reduce_mpi_cpu( ny, &h_accumulator2[0], &h_accumulator[0], comm, comm_mod, comm_mod_reduce);
    for( unsigned i=0; i<ny; i++)
        out[i] = exblas::cpu::Round( &h_accumulator[i*exblas::BIN_COUNT]);
}
#endif //MPI_VERSION

}//namespace dg