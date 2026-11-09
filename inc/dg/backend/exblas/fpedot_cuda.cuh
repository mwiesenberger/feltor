/*
 * %%%%%%%%%%%%%%%%%%%%%%%Original development%%%%%%%%%%%%%%%%%%%%%%%%%
 *  Matthias Wiesenberger, 2020, within FELTOR license
 */
/**
 *  @file fpedot_cuda.cuh
 *  @brief CUDA version of fpedot
 *
 *  @authors
 *    Developers : \n
 *        Matthias Wiesenberger -- mattwi@fysik.dtu.dk
 */
#pragma once
#include <algorithm>
#include <array>
#include "thrust/device_vector.h"
#include "accumulate.cuh"

namespace dg
{
namespace exblas{
///@cond
namespace gpu{

////////////////////////////////////////////////////////////////////////////////
// Main Kernels
////////////////////////////////////////////////////////////////////////////////
//In the first kernel each block produces exactly one FPE (because threads within a block can be synchronized and shared memory lives for a block )
//the second kernel reduces all FPEs from the first kernel (because we need a global synchronization, which is induced by separate kernel launches)

template<class T, unsigned N>
__device__ void warpReduce64( T * a, unsigned int tid, unsigned int stride)
{
    // We reduce 64 FPEs, stride is the stride to access FPEs
    // assert( 64 == 2*(max_tid+1)) i.e. tid = 0,1,....,31 i.e. a warp
    // This is a tree sum that "folds a onto itself"
    // 1. In each first iteration the 2nd half of a is summed onto the first half
    // 2. until there is only one left

    // EVERY THREAD IN A WARP SHOULD DO THE SAME INSTRUCTION (so no if(threadIdx) allowed)
    int modulo_factor = 64/2;
    while( modulo_factor >=  1 )
    {
        // Thread tid reads FPE tid+modulo_factor and writes FPE tid. For
        // tid >= modulo_factor, FPE tid is read by thread tid-modulo_factor,
        // so all threads must finish reading before any thread writes
        T y[N];
        for( unsigned k=0; k<N; k++)
            y[k] = a[k*stride+tid+modulo_factor];
        __syncwarp();
        for( unsigned k=0; k<N; k++)
        {
            T x = y[k];
            #pragma unroll
            for(unsigned i = 0; i != N; ++i) {
                T s;
                a[i*stride+tid] = KnuthTwoSum(a[i*stride+tid], x, &s);
                x = s;
            }
        }
        modulo_factor/=2;
        // We need a memory barrier between iterations
        // (because threads need to read what other threads computed previously)
        // Otherwise results are wrong!!
        // https://stackoverflow.com/questions/46467011/thread-synchronization-with-syncwarp
        __syncwarp();
    }
}

template<class T, unsigned N, unsigned THREADS_PER_BLOCK, class Functor, class ...PointerOrValues>
__global__ void fpeDOT(
    const unsigned NbElements,
    T *d_PartialFPEs,
    const unsigned stride, // of the components of d_PartialFPEs
    Functor f,
    PointerOrValues ...d_xs
) {
    // Dynamic shared memory
    // https://developer.nvidia.com/blog/using-shared-memory-cuda-cc/
    // and
    // https://stackoverflow.com/questions/27570552/templated-cuda-kernel-with-dynamic-shared-memory
    //extern __shared__ T l_fpe[]; //shared variables live for a thread block
    extern __shared__ unsigned char memory[];
    T *l_fpe = reinterpret_cast<T *>(memory);
    T *a = l_fpe + threadIdx.x;
    //Initialize FPEs
    for (unsigned i = 0; i < N; i++)
        a[i * THREADS_PER_BLOCK] = T(0);
    __syncthreads();

    //Read data from global memory and accumulate to sub-FPEs
    for(unsigned pos = blockIdx.x*THREADS_PER_BLOCK+threadIdx.x; pos < NbElements; pos += gridDim.x*THREADS_PER_BLOCK) {
        T x = f(get_element(d_xs,pos)...);

        #pragma unroll
        for(unsigned i = 0; i != N; ++i) {
            T s;
            a[i*THREADS_PER_BLOCK] = KnuthTwoSum(a[i*THREADS_PER_BLOCK], x, &s);
            x = s;
        }
    }
    __syncthreads();

    //Tree reduction ("fold FPEs on themselves")
    int modulo_factor = THREADS_PER_BLOCK/2;
    while( modulo_factor >=  64 )
    {
        if(threadIdx.x < modulo_factor) {
            for( unsigned k=0; k<N; k++)
            {
                T x = a[k*THREADS_PER_BLOCK+modulo_factor];
                #pragma unroll
                for(unsigned i = 0; i != N; ++i) {
                    T s;
                    a[i*THREADS_PER_BLOCK] = KnuthTwoSum(a[i*THREADS_PER_BLOCK], x, &s);
                    x = s;
                }
            }
        }
        modulo_factor/=2;
        __syncthreads();
    }
    // We are now left with 64 FPEs in the block
    //Now merge sub-FPEs within each warp
    if( threadIdx.x < 32) warpReduce64<T, N>( l_fpe, threadIdx.x, THREADS_PER_BLOCK);
    if( threadIdx.x == 0)
    {
        for( unsigned i=0; i<N; i++)
            d_PartialFPEs[blockIdx.x + i*stride] = l_fpe[i*THREADS_PER_BLOCK];
    }

}

////////////////////////////////////////////////////////////////////////////////
// Merging
////////////////////////////////////////////////////////////////////////////////

// ATTENTION The following values were only measured on a TITAN Xp and should
// be reevaluated on an A100. Only the cap on the number of blocks adapts to
// the device at runtime
static constexpr unsigned FPE_BLOCK_THREADS  = 64;   //# threads per block >= 64 and a power of 2
static constexpr unsigned FPE_MERGE_SIZE     = 64;   //# of FPEs a merge block reduces (fixed by warpReduce64)
static constexpr unsigned FPE_SMALL_ITEMS    = 2;    //# of elements per thread for up to FPE_MERGE_SIZE blocks
static constexpr unsigned FPE_MIN_ITEMS      = 32;   //# of elements per thread before we add blocks beyond FPE_MERGE_SIZE
static constexpr unsigned MAX_PARTIAL_FPES   = FPE_MERGE_SIZE*FPE_MERGE_SIZE; //what the two merge kernels can handle

// Every block merges up to FPE_MERGE_SIZE consecutive FPEs (of count in total)
// to one FPE. There may only be blocks with 32 threads
template<class T, unsigned N>
__global__
void fpeDOTMerge(
     T *d_fpes,
     T *d_merged,
     const unsigned count,
     const unsigned stride,        // of the components of d_fpes
     const unsigned merged_stride  // of the components of d_merged
) {
    const unsigned tid = threadIdx.x, first = blockIdx.x*FPE_MERGE_SIZE;
    T *a = d_fpes + first;
    // warpReduce64 reduces exactly 64 FPEs so we zero the unused ones
    const unsigned items = count - first;
    for( unsigned pos : {tid, tid + FPE_MERGE_SIZE/2})
        if( pos >= items)
            for( unsigned k=0; k<N; k++)
                a[k*stride + pos] = T(0);
    __syncwarp();
    warpReduce64<T, N>( a, tid, stride);
    if( tid == 0)
    {
        for( unsigned i=0; i<N; i++)
            d_merged[blockIdx.x + i*merged_stride] = a[i*stride];
    }
}

// Every block has a fixed cost and more than FPE_MERGE_SIZE blocks need an
// extra merge kernel. So use one block per FPE_SMALL_ITEMS*FPE_BLOCK_THREADS
// elements up to FPE_MERGE_SIZE blocks, then stay there until every thread has
// FPE_MIN_ITEMS elements. One merge kernel is thus enough up to
// FPE_MERGE_SIZE*FPE_BLOCK_THREADS*FPE_MIN_ITEMS = 131072 elements. More blocks
// than the device can run at the same time only wait
inline unsigned num_fpe_blocks( unsigned size, unsigned max_blocks)
{
    auto ceil_div = []( unsigned a, unsigned b){ return a/b + (a%b != 0);};
    unsigned blocks = ceil_div( size, FPE_SMALL_ITEMS*FPE_BLOCK_THREADS);
    if( blocks > FPE_MERGE_SIZE)
        blocks = std::max( FPE_MERGE_SIZE, ceil_div( size,
            FPE_BLOCK_THREADS*FPE_MIN_ITEMS));
    return std::min( std::max( blocks, 1u), max_blocks);
}

template<class T, unsigned ID>
inline T* device_fpe_buffer( unsigned size)
{
    static thrust::device_vector<T> buffer;
    if( buffer.size() < size)
        buffer.resize( size);
    return thrust::raw_pointer_cast( buffer.data());
}

// The result lives in pinned host memory (see HostResult in exdot_cuda.cuh)
template<class T, size_t N>
inline T* host_fpe()
{
    static T* fpe = []()
    {
        void* ptr = nullptr;
        cudaError_t code = cudaHostAlloc( &ptr, N*sizeof(T), cudaHostAllocDefault);
        if( code != cudaSuccess)
            throw dg::Error(dg::Message(_ping_)<<cudaGetErrorString(code));
        return static_cast<T*>(ptr);
    }();
    return fpe;
}

}//namespace gpu
///@endcond


/*!@brief GPU version of fpe generalized dot product
 *
 * @copydetails fpedot_cpu
 * @copydoc hide_gpu_static_buffer
 */
template<class T, size_t N, class Functor, class ...PointerOrValues>
__host__
void fpedot_gpu(unsigned size, std::array<T,N>& fpe, Functor f, PointerOrValues ...xs_ptr)
{
    constexpr unsigned threads = gpu::FPE_BLOCK_THREADS;
    constexpr size_t shared_memory = N*threads*sizeof(T);
    static const unsigned max_blocks = gpu::max_active_blocks(
        static_cast<void(*)(const unsigned, T*, const unsigned, Functor, PointerOrValues...)>(
            gpu::fpeDOT<T, N, threads, Functor, PointerOrValues...>),
        threads, shared_memory, gpu::MAX_PARTIAL_FPES);
    const unsigned stride = std::max( max_blocks, gpu::FPE_MERGE_SIZE);
    T *d_PartialFPEs = gpu::device_fpe_buffer<T,0>( N*stride);
    const unsigned blocks = gpu::num_fpe_blocks( size, max_blocks);
    T *h_fpe = gpu::host_fpe<T,N>();
    gpu::fpeDOT<T, N, threads, Functor, PointerOrValues...>
        <<<blocks, threads, shared_memory>>>(
            size, d_PartialFPEs, stride, f, xs_ptr...);
    if( blocks <= gpu::FPE_MERGE_SIZE)
        gpu::fpeDOTMerge<T, N><<<1, 32>>>( d_PartialFPEs, h_fpe, blocks, stride, 1);
    else
    {
        T *d_merged = gpu::device_fpe_buffer<T,1>( N*gpu::FPE_MERGE_SIZE);
        const unsigned groups = blocks/gpu::FPE_MERGE_SIZE + (blocks%gpu::FPE_MERGE_SIZE != 0);
        gpu::fpeDOTMerge<T, N><<<groups, 32>>>( d_PartialFPEs, d_merged, blocks, stride, gpu::FPE_MERGE_SIZE);
        gpu::fpeDOTMerge<T, N><<<1, 32>>>( d_merged, h_fpe, groups, gpu::FPE_MERGE_SIZE, 1);
    }
    gpu::wait_for_kernels();
    std::copy_n( h_fpe, N, fpe.begin());
}

}//namespace exblas
} //namespace dg
