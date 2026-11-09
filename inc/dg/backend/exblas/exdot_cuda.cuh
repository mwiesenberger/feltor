/*
 * %%%%%%%%%%%%%%%%%%%%%%%Original development%%%%%%%%%%%%%%%%%%%%%%%%%
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 *  All rights reserved.
 * %%%%%%%%%%%%%%%%%%%%%%%Modifications and further additions%%%%%%%%%%
 *  Matthias Wiesenberger, 2017, within FELTOR and EXBLAS licenses
 */
/**
 *  @file exdot_cuda.cuh
 *  @brief CUDA version of exdot
 *
 *  @authors
 *    Developers : \n
 *        Roman Iakymchuk  -- roman.iakymchuk@lip6.fr \n
 *        Sylvain Collange -- sylvain.collange@inria.fr \n
 *        Matthias Wiesenberger -- mattwi@fysik.dtu.dk
 */
#pragma once
#include <algorithm>
#include "thrust/device_vector.h"
#include "../exceptions.h"
#include "accumulate.cuh"

namespace dg
{
namespace exblas{
///@cond
namespace gpu{

////////////////////////////////////////////////////////////////////////////////
// Main Kernels
////////////////////////////////////////////////////////////////////////////////
//In the first kernel each block produces exactly one superacc (because threads within a block can be synchronized and shared memory lives for a block )
//the second kernel reduces all superaccs from the first kernel (because we need a global synchronization, which is induced by separate kernel launches)

template<unsigned NBFPE, unsigned WARP_COUNT, class PointerOrValue1, class PointerOrValue2>
__global__ void ExDOT(
    int64_t *d_PartialSuperaccs,
    PointerOrValue1 d_a,
    PointerOrValue2 d_b,
    const unsigned NbElements,
    volatile bool* error
) {
    __shared__ int64_t l_sa[WARP_COUNT * BIN_COUNT]; //shared variables live for a thread block (39 rows, 16 columns!)
    int64_t *l_workingBase = l_sa + (threadIdx.x & (WARP_COUNT - 1)); //the bitwise & with 15 is a modulo operation: threadIdx.x % 16
    //Initialize superaccs
    for (unsigned i = 0; i < BIN_COUNT; i++)
        l_workingBase[i * WARP_COUNT] = 0;
    __syncthreads(); //syncs all threads in a block (but not across blocks)

    //Read data from global memory and scatter it to sub-superaccs
    double a[NBFPE] = {0.0};
    for(unsigned pos = blockIdx.x*blockDim.x+threadIdx.x; pos < NbElements; pos += gridDim.x*blockDim.x) {
        //double r = 0.0;
        //double x = TwoProductFMA(get_element(d_a,pos), get_element(d_b,pos), &r);
        double x = (double)get_element(d_a,pos)*(double)get_element(d_b,pos);
        //we do not accumulate the rest of this multiplication

        //Check if the input is sane
        if( !isfinite(x) ) *error = true;

        #pragma unroll
        for(unsigned i = 0; i != NBFPE; ++i) {
            double s;
            a[i] = KnuthTwoSum(a[i], x, &s);
            x = s;
        }
        if (x != 0.0) {
            Accumulate(l_workingBase, x, WARP_COUNT);
            // Flush FPEs to superaccs
            #pragma unroll
            for(unsigned i = 0; i != NBFPE; ++i) {
                Accumulate(l_workingBase, a[i], WARP_COUNT);
                a[i] = 0.0;
            }
        }

        //if (r != 0.0) {//add the rest r in the same manner
        //    #pragma unroll
        //    for(unsigned i = 0; i != NBFPE; ++i) {
        //        double s;
        //        a[i] = KnuthTwoSum(a[i], r, &s);
        //        r = s;
        //    }
        //    if (r != 0.0) {
        //        Accumulate(l_workingBase, r, WARP_COUNT);
        //        // Flush FPEs to superaccs
        //        #pragma unroll
        //        for(unsigned i = 0; i != NBFPE; ++i) {
        //            Accumulate(l_workingBase, a[i], WARP_COUNT);
        //            a[i] = 0.0;
        //        }
        //    }
        //}
    }
    //Flush FPEs to superaccs
    #pragma unroll
    for(unsigned i = 0; i != NBFPE; ++i)
        Accumulate(l_workingBase, a[i], WARP_COUNT);
    __syncthreads();

    //Merge sub-superaccs into work-group partial-accumulator ( ATTENTION: PartialSuperacc is transposed!)
    unsigned pos = threadIdx.x;
    if(pos < WARP_COUNT) {
        int imin = IMIN, imax = IMAX;
        Normalize( l_workingBase, imin, imax, WARP_COUNT);
    }
    __syncthreads();

    if (pos < BIN_COUNT) {
        int64_t sum = 0;

        for(unsigned i = 0; i < WARP_COUNT; i++)
            sum += l_sa[pos * WARP_COUNT + i];

        d_PartialSuperaccs[blockIdx.x * BIN_COUNT + pos] = sum;
    }

    __syncthreads();
    if (pos == 0) {
        int imin = IMIN, imax = IMAX;
        Normalize(&d_PartialSuperaccs[blockIdx.x * BIN_COUNT], imin, imax);
    }
}

template<unsigned NBFPE, unsigned WARP_COUNT, class PointerOrValue1, class PointerOrValue2, class PointerOrValue3>
__global__ void ExDOT(
    int64_t *d_PartialSuperaccs,
    PointerOrValue1 d_a,
    PointerOrValue2 d_b,
    PointerOrValue3 d_c,
    const unsigned NbElements,
    volatile bool *error
) {
    __shared__ int64_t l_sa[WARP_COUNT * BIN_COUNT]; //shared variables live for a thread block (39 rows, 16 columns!)
    int64_t *l_workingBase = l_sa + (threadIdx.x & (WARP_COUNT - 1)); //the bitwise & with 15 is a modulo operation: threadIdx.x % 16
    //Initialize superaccs
    for (unsigned i = 0; i < BIN_COUNT; i++)
        l_workingBase[i * WARP_COUNT] = 0;
    __syncthreads(); //syncs all threads in a block (but not across blocks)

    //Read data from global memory and scatter it to sub-superaccs
    double a[NBFPE] = {0.0};
    for(unsigned pos = blockIdx.x*blockDim.x+threadIdx.x; pos < NbElements; pos += gridDim.x*blockDim.x) {
        //double x2 = d_a[pos]*d_c[pos]*d_b[pos];
        //double r  = 0.0, r2 = 0.0;
        //double x  = TwoProductFMA(d_a[pos], d_b[pos], &r);
        //double x2 = TwoProductFMA(x , d_c[pos], &r2);
        double x1 = (double)get_element(d_a,pos)*(double)get_element(d_b,pos);
        double x2 = x1                          *(double)get_element(d_c,pos);

        //Check if the input is sane
        if( !isfinite(x2) ) *error = true;

        if( x2 != 0.0 ) {//accumulate x2
            #pragma unroll
            for(unsigned i = 0; i != NBFPE; ++i) {
                double s;
                a[i] = KnuthTwoSum(a[i], x2, &s);
                x2 = s;
            }
            if (x2 != 0.0) {
                Accumulate(l_workingBase, x2, WARP_COUNT);
                // Flush FPEs to superaccs
                #pragma unroll
                for(unsigned i = 0; i != NBFPE; ++i) {
                    Accumulate(l_workingBase, a[i], WARP_COUNT);
                    a[i] = 0.0;
                }
            }
        }
        //if (r2 != 0.0) {//add the rest r2
        //    #pragma unroll
        //    for(unsigned i = 0; i != NBFPE; ++i) {
        //        double s;
        //        a[i] = KnuthTwoSum(a[i], r2, &s);
        //        r2 = s; //error was here r = s
        //    }
        //    if (r2 != 0.0) { //error was here r != 0.0
        //        Accumulate(l_workingBase, r2, WARP_COUNT);
        //        // Flush FPEs to superaccs
        //        #pragma unroll
        //        for(unsigned i = 0; i != NBFPE; ++i) {
        //            Accumulate(l_workingBase, a[i], WARP_COUNT);
        //            a[i] = 0.0;
        //        }
        //    }
        //}

        //if (r != 0.0) {//add the rest r*c in the same manner
        //    r2 = 0.0;
        //    x2 = TwoProductFMA(r , d_c[pos], &r2);
        //    if( x2 != 0.0) {//accumulate x2
        //        #pragma unroll
        //        for(unsigned i = 0; i != NBFPE; ++i) {
        //            double s;
        //            a[i] = KnuthTwoSum(a[i], x2, &s);
        //            x2 = s;
        //        }
        //        if (x2 != 0.0) {
        //            Accumulate(l_workingBase, x2, WARP_COUNT);
        //            // Flush FPEs to superaccs
        //            #pragma unroll
        //            for(unsigned i = 0; i != NBFPE; ++i) {
        //                Accumulate(l_workingBase, a[i], WARP_COUNT);
        //                a[i] = 0.0;
        //            }
        //        }
        //    }
        //    if (r2 != 0.0) {//add the rest r2
        //        #pragma unroll
        //        for(unsigned i = 0; i != NBFPE; ++i) {
        //            double s;
        //            a[i] = KnuthTwoSum(a[i], r2, &s);
        //            r2 = s; //error was here r = s
        //        }
        //        if (r2 != 0.0) { //error was here r != 0.0
        //            Accumulate(l_workingBase, r2, WARP_COUNT);
        //            // Flush FPEs to superaccs
        //            #pragma unroll
        //            for(unsigned i = 0; i != NBFPE; ++i) {
        //                Accumulate(l_workingBase, a[i], WARP_COUNT);
        //                a[i] = 0.0;
        //            }
        //        }
        //    }
        //}
    }
	//Flush FPEs to superaccs
    #pragma unroll
    for(unsigned i = 0; i != NBFPE; ++i)
        Accumulate(l_workingBase, a[i], WARP_COUNT);
    __syncthreads();

    //Merge sub-superaccs into work-group partial-accumulator ( ATTENTION: PartialSuperacc is transposed!)
    unsigned pos = threadIdx.x;
    if(pos<WARP_COUNT){
            int imin = IMIN, imax = IMAX;
            Normalize(l_workingBase, imin, imax, WARP_COUNT);
        }
    __syncthreads();
    if (pos < BIN_COUNT) {
        int64_t sum = 0;

        for(unsigned i = 0; i < WARP_COUNT; i++)
            sum += l_sa[pos * WARP_COUNT + i];

        d_PartialSuperaccs[blockIdx.x * BIN_COUNT + pos] = sum;
    }

    __syncthreads();
    if (pos == 0) {
        int imin = IMIN, imax = IMAX;
        Normalize(&d_PartialSuperaccs[blockIdx.x * BIN_COUNT], imin, imax);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Merging
////////////////////////////////////////////////////////////////////////////////
////////////// parameters for Kernel execution            //////////////////////
static constexpr unsigned WARP_COUNT            = 16;  //# of sub superaccs per block in ExDOT
static constexpr unsigned EXDOT_BLOCK_THREADS   = 128; //# threads per block in ExDOT
static constexpr unsigned MIN_ITEMS_PER_THREAD  = 16;  //# of elements per thread before we add blocks beyond MERGE_SIZE
static constexpr unsigned MERGE_SIZE            = 64;  //# of superaccs a block sums; <= 128 since normalized words are < 2^56
static constexpr unsigned MERGE_BLOCK_THREADS   = 64;  //we need minimum 39 of those
static constexpr unsigned MAX_PARTIAL_SUPERACCS = MERGE_SIZE*MERGE_SIZE; //what the two merge kernels can handle

template<unsigned MERGE>
__global__
void ExDOTComplete(
     int64_t *d_PartialSuperaccs,
     const unsigned count
) {
    unsigned lid = threadIdx.x;
    unsigned gid = blockIdx.x;

    if (lid < BIN_COUNT) { //every block sums its assigned superaccs
        int64_t sum = 0;

        for(unsigned i = 0; i < MERGE && gid * MERGE + i < count; i++)
            sum += d_PartialSuperaccs[(gid * MERGE + i) * BIN_COUNT + lid];

        d_PartialSuperaccs[gid * BIN_COUNT * MERGE + lid] = sum;
    }

    __syncthreads(); //syncs all threads in a block (but not across blocks)
    if (lid == 0) { //every block normalize its summed superacc
        int imin = IMIN, imax = IMAX;
        Normalize(&d_PartialSuperaccs[gid * BIN_COUNT * MERGE], imin, imax);
    }
}
//MW: global synchronization here through separate Kernels
//one block of threads with at least 39 threads
//The result is not normalized: every user of a superaccumulator has to do that
template<unsigned STRIDE>
__global__
void ExDOTCompleteFinal(
     int64_t *d_PartialSuperaccs,
     int64_t *d_superacc,
     const unsigned count // at most MERGE_SIZE
) {
    unsigned lid = threadIdx.x;
    if (lid < BIN_COUNT) {
        int64_t sum = 0;

        for(unsigned i = 0; i < count; i++)
            sum += d_PartialSuperaccs[i * BIN_COUNT * STRIDE + lid];

        d_superacc[lid] = sum;
    }
}

// How many blocks of ExDOT the device can run at the same time. Queried from
// the runtime, so that the number of blocks adapts to the device in use
template<class Kernel>
inline unsigned max_active_blocks( Kernel kernel)
{
    int device = 0, blocks_per_sm = 0;
    cudaDeviceProp prop;
    cudaError_t code = cudaGetDevice( &device);
    if( code == cudaSuccess)
        code = cudaGetDeviceProperties( &prop, device);
    if( code == cudaSuccess)
        code = cudaOccupancyMaxActiveBlocksPerMultiprocessor( &blocks_per_sm,
            kernel, EXDOT_BLOCK_THREADS, 0);
    if( code != cudaSuccess)
        throw dg::Error(dg::Message(_ping_)<<cudaGetErrorString(code));
    unsigned blocks = std::max( prop.multiProcessorCount*blocks_per_sm, 1);
    return std::min( blocks, MAX_PARTIAL_SUPERACCS);
}

// Every block has a fixed cost and more than MERGE_SIZE blocks need an extra
// merge kernel. So use one block per EXDOT_BLOCK_THREADS elements up to
// MERGE_SIZE blocks, then stay there until every thread has MIN_ITEMS_PER_THREAD
// elements. More blocks than the device can run at the same time only wait
inline unsigned num_exdot_blocks( unsigned size, unsigned max_blocks)
{
    auto ceil_div = []( unsigned a, unsigned b){ return a/b + (a%b != 0);};
    unsigned blocks = ceil_div( size, EXDOT_BLOCK_THREADS);
    if( blocks > MERGE_SIZE)
        blocks = std::max( MERGE_SIZE, ceil_div( size,
            EXDOT_BLOCK_THREADS*MIN_ITEMS_PER_THREAD));
    return std::min( std::max( blocks, 1u), max_blocks);
}

inline int64_t* partial_superaccs( unsigned count)
{
    static thrust::device_vector<int64_t> partial;
    if( partial.size() < count*BIN_COUNT)
        partial.resize( count*BIN_COUNT);
    return thrust::raw_pointer_cast( partial.data());
}

// With at most MERGE_SIZE blocks one kernel is enough, which saves a kernel
// launch. With num_exdot_blocks that is the case up to
// MERGE_SIZE*EXDOT_BLOCK_THREADS*MIN_ITEMS_PER_THREAD = 131072 elements
inline void merge_partial_superaccs( int64_t* d_partial, unsigned blocks, int64_t* superacc)
{
    if( blocks <= MERGE_SIZE)
    {
        ExDOTCompleteFinal<1><<<1, MERGE_BLOCK_THREADS>>>( d_partial, superacc, blocks);
        return;
    }
    unsigned groups = blocks/MERGE_SIZE + (blocks%MERGE_SIZE != 0);
    ExDOTComplete<MERGE_SIZE><<<groups, MERGE_BLOCK_THREADS>>>( d_partial, blocks);
    ExDOTCompleteFinal<MERGE_SIZE><<<1, MERGE_BLOCK_THREADS>>>( d_partial, superacc, groups);
}

// The result of a dot product. It lives in pinned host memory, so that the
// last kernel writes it directly to the host and no cudaMemcpy is needed.
// A single write to pinned memory is slower than to device memory, but the
// result is only BIN_COUNT words, while a cudaMemcpy has a fixed cost however
// little it moves. On a TITAN Xp the last kernel plus cudaMemcpy took 13.2us
// and the last kernel writing here plus cudaDeviceSynchronize 4.7us, which
// matters since a dot product at small sizes is dominated by such fixed costs.
// Only do this for small data that kernels write and never read: reading
// pinned memory from a kernel or writing large vectors to it is slow.
struct HostResult
{
    int64_t superacc[BIN_COUNT];
    volatile bool error; // only written by the kernels for NaN or Inf
};

inline HostResult* host_result()
{
    static HostResult* result = []()
    {
        void* ptr = nullptr;
        cudaError_t code = cudaHostAlloc( &ptr, sizeof(HostResult), cudaHostAllocDefault);
        if( code != cudaSuccess)
            throw dg::Error(dg::Message(_ping_)<<cudaGetErrorString(code));
        return static_cast<HostResult*>(ptr);
    }();
    return result;
}

// Wait for the launched kernels and copy their result
inline void finish_exdot( const HostResult* result, int64_t* h_superacc, int* status)
{
    cudaError_t code = cudaDeviceSynchronize();
    if( code != cudaSuccess)
        throw dg::Error(dg::Message(_ping_)<<cudaGetErrorString(code));
    // This test checks for errors in the current stream, the error may come
    // from any kernel prior to this point not necessarily the ones of the dot
    code = cudaGetLastError( );
    if( code != cudaSuccess)
        throw dg::Error(dg::Message(_ping_)<<cudaGetErrorString(code));
    std::copy_n( result->superacc, BIN_COUNT, h_superacc);
    *status = result->error ? 1 : 0;
}

}//namespace gpu
///@endcond

///@brief GPU version of exact dot product
///@copydoc hide_exdot2
///@copydoc hide_hostacc
///@copydoc hide_gpu_static_buffer
template<class PointerOrValue1, class PointerOrValue2, size_t NBFPE=3>
__host__
void exdot_gpu(unsigned size, PointerOrValue1 x1_ptr, PointerOrValue2 x2_ptr, int64_t* h_superacc, int* status)
{
    static_assert( has_floating_value<PointerOrValue1>::value, "PointerOrValue1 needs to be T or T* with T one of (const) float or (const) double");
    static_assert( has_floating_value<PointerOrValue2>::value, "PointerOrValue2 needs to be T or T* with T one of (const) float or (const) double");
    static const unsigned max_blocks = gpu::max_active_blocks( static_cast<void(*)(int64_t*, PointerOrValue1, PointerOrValue2, unsigned, volatile bool*)>(
        gpu::ExDOT<NBFPE, gpu::WARP_COUNT, PointerOrValue1, PointerOrValue2>));
    int64_t *d_PartialSuperaccs = gpu::partial_superaccs( max_blocks);
    const unsigned blocks = gpu::num_exdot_blocks( size, max_blocks);
    gpu::HostResult *result = gpu::host_result();
    result->error = false;
    gpu::ExDOT<NBFPE, gpu::WARP_COUNT><<<blocks, gpu::EXDOT_BLOCK_THREADS>>>( d_PartialSuperaccs, x1_ptr, x2_ptr, size, &result->error);
    gpu::merge_partial_superaccs( d_PartialSuperaccs, blocks, result->superacc);
    gpu::finish_exdot( result, h_superacc, status);
}

///@brief GPU version of exact dot product
///@copydoc hide_exdot3
///@copydoc hide_hostacc
///@copydoc hide_gpu_static_buffer
template<class PointerOrValue1, class PointerOrValue2, class PointerOrValue3, size_t NBFPE=3>
__host__
void exdot_gpu(unsigned size, PointerOrValue1 x1_ptr, PointerOrValue2 x2_ptr, PointerOrValue3 x3_ptr, int64_t* h_superacc, int* status)
{
    static_assert( has_floating_value<PointerOrValue1>::value, "PointerOrValue1 needs to be T or T* with T one of (const) float or (const) double");
    static_assert( has_floating_value<PointerOrValue2>::value, "PointerOrValue2 needs to be T or T* with T one of (const) float or (const) double");
    static_assert( has_floating_value<PointerOrValue3>::value, "PointerOrValue3 needs to be T or T* with T one of (const) float or (const) double");
    static const unsigned max_blocks = gpu::max_active_blocks( static_cast<void(*)(int64_t*, PointerOrValue1, PointerOrValue2, PointerOrValue3, unsigned, volatile bool*)>(
        gpu::ExDOT<NBFPE, gpu::WARP_COUNT, PointerOrValue1, PointerOrValue2, PointerOrValue3>));
    int64_t *d_PartialSuperaccs = gpu::partial_superaccs( max_blocks);
    const unsigned blocks = gpu::num_exdot_blocks( size, max_blocks);
    gpu::HostResult *result = gpu::host_result();
    result->error = false;
    gpu::ExDOT<NBFPE, gpu::WARP_COUNT><<<blocks, gpu::EXDOT_BLOCK_THREADS>>>( d_PartialSuperaccs, x1_ptr, x2_ptr, x3_ptr, size, &result->error);
    gpu::merge_partial_superaccs( d_PartialSuperaccs, blocks, result->superacc);
    gpu::finish_exdot( result, h_superacc, status);
}

}//namespace exblas
} //namespace dg
