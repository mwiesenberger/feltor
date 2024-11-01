#ifndef _DG_BLAS_OMP_
#define _DG_BLAS_OMP_
#include <omp.h>
#include <thrust/transform_reduce.h>
#include <thrust/system/omp/execution_policy.h>
#include "config.h"
#include "blas1_serial.h"
#include "exblas/exdot_omp.h"
namespace dg
{
namespace blas1
{
namespace detail
{
constexpr int MIN_SIZE=100;//don't parallelize if work is too small

template<class PointerOrValue1, class PointerOrValue2>
inline std::vector<int64_t> doDot_dispatch( OmpTag, unsigned size, PointerOrValue1 x_ptr, PointerOrValue2 y_ptr) {
    std::vector<int64_t> h_superacc(exblas::BIN_COUNT);
    int status = 0;
    if(size<MIN_SIZE)
        exblas::exdot_cpu( size, x_ptr,y_ptr, &h_superacc[0], &status);
    else
        exblas::exdot_omp( size, x_ptr,y_ptr, &h_superacc[0], &status);
    if(status != 0)
        throw dg::Error(dg::Message(_ping_)<<"OMP Dot failed since one of the inputs contains NaN or Inf");
    return h_superacc;
}
template<class PointerOrValue1, class PointerOrValue2, class PointerOrValue3>
inline std::vector<int64_t> doDot_dispatch( OmpTag, unsigned size, PointerOrValue1 x_ptr, PointerOrValue2 y_ptr, PointerOrValue3 z_ptr) {
    std::vector<int64_t> h_superacc(exblas::BIN_COUNT);
    int status = 0;
    if(size<MIN_SIZE)
        exblas::exdot_cpu( size, x_ptr,y_ptr,z_ptr, &h_superacc[0], &status);
    else
        exblas::exdot_omp( size, x_ptr,y_ptr,z_ptr, &h_superacc[0], &status);
    if(status != 0)
        throw dg::Error(dg::Message(_ping_)<<"OMP Dot failed since one of the inputs contains NaN or Inf");
    return h_superacc;
}

template< class Subroutine, class PointerOrValue, class ...PointerOrValues>
inline void doSubroutine_omp( int size, Subroutine f, PointerOrValue x, PointerOrValues... xs)
{
#pragma omp for nowait
    for( int i=0; i<size; i++)
        //f(x[i], xs[i]...);
        //f(thrust::raw_reference_cast(*(x+i)), thrust::raw_reference_cast(*(xs+i))...);
        f(get_element(x,i), get_element(xs,i)...);
}

template< class Subroutine, class PointerOrValue, class ...PointerOrValues>
inline void doSubroutine_dispatch( OmpTag, int size, Subroutine f, PointerOrValue x, PointerOrValues... xs)
{
    if(omp_in_parallel())
    {
        doSubroutine_omp( size, f, x, xs... );
        return;
    }
    if(size>MIN_SIZE)
    {
        #pragma omp parallel
        {
            doSubroutine_omp( size, f, x, xs...);
        }
    }
    else
        doSubroutine_dispatch( SerialTag(), size, f, x, xs...);
}

template<class T, class Pointer, class BinaryOp, class UnaryOp>
inline T doReduce_dispatch( OmpTag, int size, Pointer x, T init, BinaryOp op,
        UnaryOp unary_op)
{
    return thrust::transform_reduce(thrust::omp::par, x, x+size, unary_op, init, op);
}
template<class F, class G, size_t N, class Pointer, class ...PointerOrValues>
void doKronecker_omp( Pointer y, size_t size, F f, G g, const std::array<size_t, N>& sizes, PointerOrValues ...xs)
{
#pragma omp for nowait
    for( unsigned u=0; u<size; u++)
    {
        std::array<size_t, N> current;
        current[0] = u%sizes[0];
        size_t remain = u/sizes[0];
        for( unsigned k=1; k<N; k++)
        {
            current[k] = remain%sizes[k];
            remain = remain/sizes[k];
        }
        call_host_F( f, g, y, u, &current[0], std::make_index_sequence<N>(), xs ...);
    }
}
template<class F, class G, size_t N, class Pointer, class ...PointerOrValues>
void doKronecker_dispatch( OmpTag, Pointer y, size_t size, F f, G g, const std::array<size_t, N>& sizes, PointerOrValues ...xs)
{
    if(omp_in_parallel())
    {
        doKronecker_omp( y, size, f, g, sizes, xs... );
        return;
    }
    if(size>MIN_SIZE)
    {
        #pragma omp parallel
        {
            doKronecker_omp( y, size, f, g, sizes, xs... );
        }
    }
    else
        doKronecker_dispatch( SerialTag(), y, size, f, g, sizes, xs...);
}

}//namespace detail
}//namespace blas1
}//namespace dg
#endif //_DG_BLAS_OMP_
