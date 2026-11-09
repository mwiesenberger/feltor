#ifndef _DG_BLAS_CUDA_
#define _DG_BLAS_CUDA_
#include <thrust/transform_reduce.h>
#include <thrust/system/cuda/execution_policy.h>
#include <type_traits>
#include "exceptions.h"
#include "exblas/exdot_cuda.cuh"
#include "exblas/fpedot_cuda.cuh"
namespace dg
{
struct equals; // defined in subroutines.h
namespace blas1
{
namespace detail
{
template<class T, size_t N, class Functor, class ...PointerOrValues>
inline void doDot_fpe_dispatch( CudaTag, unsigned size, std::array<T,N>& fpe,
    Functor f, PointerOrValues ...xs_ptr)
{
    exblas::fpedot_gpu<T,N,Functor,PointerOrValues...>( size, fpe, f, xs_ptr...);
}

template<class PointerOrValue1, class PointerOrValue2>
inline std::vector<int64_t> doDot_dispatch( CudaTag, int* status, unsigned
    size, PointerOrValue1 x_ptr, PointerOrValue2 y_ptr)
{
    std::vector<int64_t> h_superacc(exblas::BIN_COUNT);
    exblas::exdot_gpu( size, x_ptr,y_ptr, &h_superacc[0], status );
    return h_superacc;
}
template<class PointerOrValue1, class PointerOrValue2, class PointerOrValue3>
inline std::vector<int64_t> doDot_dispatch( CudaTag, int* status, unsigned size, PointerOrValue1 x_ptr, PointerOrValue2 y_ptr, PointerOrValue3 z_ptr) {
    std::vector<int64_t> h_superacc(exblas::BIN_COUNT);
    exblas::exdot_gpu( size, x_ptr,y_ptr,z_ptr, &h_superacc[0], status);
    return h_superacc;
}

template<class T>
__device__
inline T get_device_element( T x, int i){
	return x;
}
template<class T>
__device__
inline T& get_device_element( T* x, int i){
	return *(x+i);
}

template<class Subroutine, class PointerOrValue, class ...PointerOrValues>
 __global__ void subroutine_kernel( int size, Subroutine f, PointerOrValue x, PointerOrValues... xs)
{
    const int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    const int grid_size = gridDim.x*blockDim.x;
    //every thread takes num_is/grid_size is
    for( int i = thread_id; i<size; i += grid_size)
        f(get_device_element(x,i), get_device_element(xs,i)...);
        //f(x[i], xs[i]...);
        //f(thrust::raw_reference_cast(*(x+i)), thrust::raw_reference_cast(*(xs+i))...);
}

template< class Subroutine, class PointerOrValue, class ...PointerOrValues>
inline void doSubroutine_dispatch( CudaTag, int size, Subroutine f, PointerOrValue x, PointerOrValues... xs)
{
    const size_t BLOCK_SIZE = 256;
    const size_t NUM_BLOCKS = std::min<size_t>((size-1)/BLOCK_SIZE+1, 65000);
    subroutine_kernel<Subroutine, PointerOrValue, PointerOrValues...><<<NUM_BLOCKS, BLOCK_SIZE>>>(size, f, x, xs...);
}

// A device to device copy of equal types needs no kernel: on a TITAN Xp
// cudaMemcpy is never slower and up to 1.9x faster for a few 1e5 elements
template< class T>
inline void doSubroutine_dispatch( CudaTag, int size, const dg::equals&, const T* x, T* y)
{
    if( x == y) // cudaMemcpy does not allow overlapping ranges
        return;
    cudaMemcpy( y, x, size_t(size)*sizeof(T), cudaMemcpyDeviceToDevice);
}

template<class T, class Pointer, class BinaryOp, class UnaryOp>
inline T doReduce_dispatch( CudaTag, int size, Pointer x, T init, BinaryOp op,
        UnaryOp unary_op)
{
    return thrust::transform_reduce(thrust::cuda::par, x, x+size, unary_op,
            init, op);
}

// Note: Here the universal reference is really important vs "F f" else we copy f on every call
template<class Binary, class F, class Pointer, std::size_t ...I, class ...PointerOrValues>
__device__
inline void call_device_F( Binary && binary, F && f, Pointer y, int i, const unsigned* a,
        std::index_sequence<I...>, PointerOrValues ... xs)
{
    binary( f( get_device_element( xs, a[I])...), y[i]);
}

// The sizes travel as a kernel argument: a device vector would cost an
// allocation per call and a memory load per element and dimension
template<size_t N>
struct KroneckerSizes
{
    unsigned s[N];
};

template<class Binary, class F, size_t N, class Pointer, class ...PointerOrValues>
__global__ void kronecker_kernel( int size, KroneckerSizes<N> sizes, Pointer y,
        Binary binary, F f, PointerOrValues ...xs)
{
    const int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    const int grid_size = gridDim.x*blockDim.x;
    unsigned current[N];
    for( int i = thread_id; i<size; i += grid_size)
    {
        // The total size is an unsigned, so 32 bit divisions suffice
        unsigned remain = i;
        for( unsigned k=0; k<N; k++)
        {
            unsigned quotient = remain/sizes.s[k];
            current[k] = remain - quotient*sizes.s[k];
            remain = quotient;
        }
        call_device_F( binary, f, y, i, current, std::make_index_sequence<N>(), xs...);
    }
}

template<class Binary, class F, size_t N, class Pointer, class ...PointerOrValues>
inline void doKronecker_dispatch( dg::CudaTag, Pointer y, size_t size, Binary &&
        binary, F && f, const std::array<size_t, N>& sizes, PointerOrValues ...xs)
{
    const size_t BLOCK_SIZE = 256;
    const size_t NUM_BLOCKS = std::min<size_t>((size-1)/BLOCK_SIZE+1, 65000);
    KroneckerSizes<N> ksizes;
    for( unsigned k=0; k<N; k++)
        ksizes.s[k] = sizes[k];
    // The kernel must take binary and f by value: Binary and F may be
    // deduced as lvalue references to host memory
    kronecker_kernel<std::decay_t<Binary>,std::decay_t<F>,N,Pointer,PointerOrValues...>
        <<<NUM_BLOCKS, BLOCK_SIZE>>>(size,
            ksizes,
            y,
            std::forward<Binary>(binary), std::forward<F>(f),
            xs...);
}


}//namespace detail
}//namespace blas1
}//namespace dg
#endif //_DG_BLAS_CUDA_
