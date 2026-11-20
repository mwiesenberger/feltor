#ifndef _DG_BLAS_OMP_
#define _DG_BLAS_OMP_
#include <omp.h>
#include <thrust/transform_reduce.h>
#include <thrust/system/omp/execution_policy.h>
#include "config.h"
#include "blas1_serial.h"
#include "exblas/exdot_omp.h"
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
namespace dg
{
namespace blas1
{
namespace detail
{
// Don't parallelize if work is too small: a parallel region costs about 5us
// on a 6 core i5-14400F. 512 slowed none of the measured kernels there, while
// 4096 made expensive functors and stencils slower.
constexpr int MIN_SIZE=512;

// The sum of the integral values f( xs[i]...); integer addition does not
// depend on the order
template<class Functor, class ...PointerOrValues>
inline auto doDot_integral_dispatch( OmpTag, unsigned size, Functor f,
    PointerOrValues ...xs)
{
    if( size < MIN_SIZE)
        return doDot_integral_dispatch( SerialTag(), size, f, xs...);
    using T = exblas::get_result_type<Functor, PointerOrValues...>;
    return thrust::transform_reduce( thrust::omp::par,
        thrust::counting_iterator<unsigned>(0),
        thrust::counting_iterator<unsigned>(size),
        [=]( unsigned i) mutable { return f( get_element( xs, i)...);},
        T(0), thrust::plus<T>());
}

// The superaccumulators of the exact sum of f( xs[i]...)
template<class Functor, class ...PointerOrValues>
inline std::vector<int64_t> doDot_dispatch( OmpTag, int * status, unsigned size,
    Functor f, PointerOrValues ...xs)
{
    constexpr unsigned K = exblas::ValueTraits<exblas::get_result_type<Functor, PointerOrValues...>>::K;
    std::vector<int64_t> h_superacc(K*exblas::BIN_COUNT);
    exblas::exsum_omp( size, &h_superacc[0], status, f, xs...);
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
// Called by every thread of a parallel region. The serial kronecker_range in
// blas1_serial.h computes any interval of the result, starting and ending in
// the middle of a row if needed (see the design comment there). Each thread
// therefore takes a contiguous share of size/tnum elements and runs the
// serial, vectorised row loop on it. f and g are taken by value, so each
// thread works on its own copy of the functors.
template<class F, class G, size_t N, class Pointer, class ...PointerOrValues>
void doKronecker_omp( Pointer y, size_t size, F f, G g, const std::array<size_t, N>& sizes, PointerOrValues ...xs)
{
    // Split by elements, not rows, so that few long rows still use all threads
    size_t tid = omp_get_thread_num();
    size_t tnum = omp_get_num_threads();
    kronecker_range( y, tid*size/tnum, (tid+1)*size/tnum, f, g, sizes, xs...);
}
template<class F, class G, size_t N, class Pointer, class ...PointerOrValues>
void doKronecker_dispatch( OmpTag, Pointer y, size_t size, F&& f, G&& g, const std::array<size_t, N>& sizes, PointerOrValues ...xs)
{
    // f and g are copied, not forwarded: every thread calls doKronecker_omp,
    // and a forwarded rvalue would be moved from by more than one thread
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
        doKronecker_dispatch( SerialTag(), y, size, std::forward<F>(f), std::forward<G>(g), sizes, xs...);
}

}//namespace detail
}//namespace blas1
}//namespace dg
#endif //_DG_BLAS_OMP_
