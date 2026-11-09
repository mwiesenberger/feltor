#ifndef _DG_BLAS_SERIAL_
#define _DG_BLAS_SERIAL_
#include <algorithm>
#include <cstring>
#include "config.h"
#include "exceptions.h"
#include "execution_policy.h"
#include "exblas/exdot_serial.h"
#include "exblas/fpedot_serial.h"

namespace dg
{
struct equals; // defined in subroutines.h
namespace blas1
{
namespace detail
{
template<class T, size_t N, class Functor, class ...PointerOrValues>
inline void doDot_fpe_dispatch( SerialTag, unsigned size, std::array<T,N>& fpe,
    Functor f, PointerOrValues ...xs_ptr)
{
    exblas::fpedot_cpu<T,N,Functor,PointerOrValues...>( size, fpe, f, xs_ptr...);
}
template<class PointerOrValue1, class PointerOrValue2>
inline std::vector<int64_t> doDot_dispatch( SerialTag, int* status, unsigned size,
    PointerOrValue1 x_ptr, PointerOrValue2 y_ptr)
{
    std::vector<int64_t> h_superacc(exblas::BIN_COUNT);
    exblas::exdot_cpu( size, x_ptr,y_ptr, &h_superacc[0], status) ;
    return h_superacc;
}
template<class PointerOrValue1, class PointerOrValue2, class PointerOrValue3>
inline std::vector<int64_t> doDot_dispatch( SerialTag, int* status, unsigned size,
    PointerOrValue1 x_ptr, PointerOrValue2 y_ptr, PointerOrValue3 z_ptr)
{
    std::vector<int64_t> h_superacc(exblas::BIN_COUNT);
    exblas::exdot_cpu( size, x_ptr,y_ptr,z_ptr, &h_superacc[0], status) ;
    return h_superacc;
}

template<class T>
inline T get_element( T x, int){
	return x;
}
template<class T>
inline T& get_element( T* x, int i){
	return *(x+i);
}
template< class Subroutine, class PointerOrValue, class ...PointerOrValues>
inline void doSubroutine_dispatch( SerialTag, int size, Subroutine f, PointerOrValue x, PointerOrValues... xs)
{
    for( int i=0; i<size; i++)
    {
        f(get_element(x,i), get_element(xs,i)...);
        //f(x[i], xs[i]...);
        //f(thrust::raw_reference_cast(*(x+i)), thrust::raw_reference_cast(*(xs+i))...);
    }
}

// A copy of equal types is a memcpy: on an i5-14400F it is up to 2x faster
// than the loop for vectors larger than the caches
template< class T>
inline void doSubroutine_dispatch( SerialTag, int size, const dg::equals&, const T* x, T* y)
{
    if( x == y) // memcpy does not allow overlapping ranges
        return;
    std::memcpy( y, x, size_t(size)*sizeof(T));
}

template<class T, class Pointer, class BinaryOp, class UnaryOp>
inline T doReduce_dispatch( SerialTag, int size, Pointer x, T init, BinaryOp
        op, UnaryOp unary_op)
{
    for(int i=0; i<size; i++)
        init = op( init, unary_op(x[i]));
    return init;
}
// Design of the CPU kronecker product
// y_u = f( x0_{i0}, x1_{i1}, ..., x_{N-1, i_{N-1}}) with the index
// u = i0 + n0*(i1 + n1*(i2 + ...)), i.e. x0 varies fastest.
// The result therefore consists of "rows" of n0 consecutive elements in which
// only i0 changes while i1, ..., i_{N-1} are fixed. The implementation has
// three layers:
// 1. kronecker_range computes y_u for u in an arbitrary interval [begin, end).
//    It computes the multi-index of begin once with divisions and then walks
//    row by row, advancing i1, i2, ... with a counter (an odometer) instead of
//    divisions per element.
// 2. The kronecker_row overload with the index_sequence reads the elements
//    x1_{i1}, ..., x_{N-1, i_{N-1}} of the current row once, so that they
//    become loop invariant values instead of loads inside the loop.
// 3. The innermost kronecker_row is a plain loop over i0 with unit stride in
//    x0 and y, which the compiler vectorises like the subroutine loop.
// Because the interval is arbitrary, a range may start and end in the middle
// of a row. The OpenMP backend (blas1_omp.h) uses this: each thread calls
// kronecker_range on its own share of elements, so the work is evenly split
// even if there are fewer rows than threads (e.g. a single long row).
// Scalars are allowed among the xs: get_element ignores the index for them.

// Note: Here the universal reference is really important vs "F f" else we copy f on every call
// y_{u+j} = f( x0_{i0+j}, xs...) for j < length, where xs are the elements
// of the other arguments in this row; a plain loop the compiler can vectorise
template<class B, class F, class Pointer, class PointerOrValue0, class ...Values>
inline void kronecker_row( B&& binary, F&& f, Pointer y, size_t u, size_t i0,
        size_t length, PointerOrValue0 x0, Values ... xs)
{
    // int like get_element, else the compiler cannot see the unit stride
    for( int j=0; j<(int)length; j++)
        binary( f( get_element( x0, (int)i0+j), xs...), y[u+j]);
}
// Reads the row's elements of x1, ..., x_{N-1} at the indices current[1],
// ..., current[N-1] (I runs over 0, ..., N-2) and passes them on by value
template<class B, class F, class Pointer, size_t ...I, class PointerOrValue0, class ...PointerOrValues>
inline void kronecker_row( B&& binary, F&& f, Pointer y, size_t u, size_t length,
        const size_t* current, std::index_sequence<I...>, PointerOrValue0 x0,
        PointerOrValues ...xs)
{
    kronecker_row( binary, f, y, u, current[0], length, x0,
            get_element( xs, current[I+1])...);
}
// y_u = f( x0_{i0}, x1_{i1}, ...) for u in [begin, end), where
// u = i0 + n0*(i1 + n1*(i2 + ...)): rows of n0 elements, in which only i0 changes
// binary and f are passed on by name, not forwarded, since they are used once
// per row
template<class B, class F, size_t N, class Pointer, class ...PointerOrValues>
void kronecker_range( Pointer y, size_t begin, size_t end, B&& binary, F&& f,
        const std::array<size_t, N>& sizes, PointerOrValues ...xs)
{
    if( begin >= end)
        return;
    // Multi-index (i0, i1, ...) of begin; the only divisions in the function
    std::array<size_t, N> current;
    size_t remain = begin;
    for( unsigned k=0; k<N; k++)
    {
        current[k] = remain%sizes[k];
        remain = remain/sizes[k];
    }
    for( size_t u = begin; u < end; )
    {
        // The first and the last row of the range may be partial
        size_t length = std::min( end - u, sizes[0] - current[0]);
        kronecker_row( binary, f, y, u, length, &current[0],
            std::make_index_sequence<N-1>(), xs...);
        u += length;
        // Next row: i0 starts at 0 and (i1, i2, ...) advance like an odometer
        current[0] = 0;
        for( unsigned k=1; k<N; k++)
        {
            current[k] ++;
            if( current[k] == sizes[k])
                current[k] = 0;
            else
                break;
        }
    }
}
template<class B, class F, size_t N, class Pointer, class ...PointerOrValues>
void doKronecker_dispatch( dg::SerialTag, Pointer y, size_t size, B&& binary, F&& f, const std::array<size_t, N>& sizes, PointerOrValues ...xs)
{
    kronecker_range( y, 0, size, binary, f, sizes, xs...);
}


}//namespace detail
}//namespace blas1
}//namespace dg
#endif //_DG_BLAS_SERIAL_
