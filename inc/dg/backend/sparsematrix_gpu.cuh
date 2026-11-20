#pragma once

#include <vector>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "config.h"

namespace dg
{
namespace detail
{
// The GPU matrix-vector product sums each row in exactly the order of
// spmv_cpu_kernel, so the result is bitwise equal to the CPU and OpenMP
// versions. The sum of one row is a sequential chain of FMAs and cannot be
// split across threads. A thread block therefore first loads alpha*A_ij and
// x_j of a group of consecutive rows coalesced into shared memory, then one
// thread per row runs the chain. A row with more than CSR_GPU_NNZ_PER_BLOCK
// entries gets a block of its own and is summed in chunks of that size, again
// in the order of the CPU. As on the CPU y is not read for beta == 0.
//
// The row blocks and the coalesced loads are CSR-stream from J. L. Greathouse
// and M. Daga, "Efficient Sparse Matrix-Vector Multiplication on GPUs using
// the CSR Storage Format", SC14 (2014), doi:10.1109/SC.2014.68. Unlike there,
// the rows are not reduced in parallel and long rows do not use CSR-vector,
// because both would change the order of summation
//
// The two constants were measured on a TITAN Xp only (FCI interpolation,
// projection and prolongation matrices, float and double): the load loop must
// have a fixed trip count CSR_GPU_NNZ_PER_BLOCK / CSR_GPU_BLOCK_SIZE so that
// all loads of a thread are in flight at once
static constexpr int CSR_GPU_BLOCK_SIZE = 128;
static constexpr int CSR_GPU_NNZ_PER_BLOCK = 1024;

// Partition of the rows into thread blocks: block b computes rows
// row_blocks[b] to row_blocks[b+1]. It depends only on the row offsets, i.e.
// the sparsity pattern
struct CSRCache_gpu
{
    void forget()
    {
        m_active = false;
        m_row_blocks.clear();
    }
    bool isUpToDate() const { return m_active;}
    // pos is the device pointer to the num_rows + 1 row offsets
    template<class I>
    void update( size_t num_rows, const I* pos)
    {
        thrust::host_vector<I> h_pos( thrust::device_pointer_cast( pos),
            thrust::device_pointer_cast( pos) + num_rows + 1);
        std::vector<int> row_blocks{0};
        int r = 0;
        while( r < (int)num_rows)
        {
            const int r0 = r;
            if( h_pos[r0+1] - h_pos[r0] > CSR_GPU_NNZ_PER_BLOCK)
                r = r0+1;
            else
                while( r < (int)num_rows && h_pos[r+1] - h_pos[r0] <= CSR_GPU_NNZ_PER_BLOCK)
                    r++;
            row_blocks.push_back( r);
        }
        m_row_blocks.assign( row_blocks.begin(), row_blocks.end());
        m_active = true;
    }
    const int* row_blocks() const { return thrust::raw_pointer_cast( m_row_blocks.data());}
    size_t num_blocks() const { return m_row_blocks.size() - 1;}

    private:
    bool m_active = false;
    thrust::device_vector<int> m_row_blocks;
};

template<class I, class V, class value_type, class C1, class C2>
__global__ void spmv_csr_stream_kernel( const int* RESTRICT row_blocks,
    const I* RESTRICT A_pos, const I* RESTRICT A_idx, const V* RESTRICT A_val,
    value_type alpha, value_type beta, const C1* RESTRICT x_ptr, C2* RESTRICT y_ptr)
{
    constexpr int NB = CSR_GPU_NNZ_PER_BLOCK, BS = CSR_GPU_BLOCK_SIZE;
    __shared__ value_type s_ax[NB]; // alpha*A_ij
    __shared__ value_type s_x[NB];  // x_j
    const int r0 = row_blocks[blockIdx.x], r1 = row_blocks[blockIdx.x+1];
    const I base = A_pos[r0], end = A_pos[r1];
    if( end - base <= NB)
    {
        const int len = end - base;
#pragma unroll
        for( int u = 0; u < NB/BS; u++)
        {
            const int k = threadIdx.x + u*BS;
            if( k < len)
            {
                s_ax[k] = alpha*A_val[base+k];
                s_x[k] = x_ptr[A_idx[base+k]];
            }
        }
        __syncthreads();
        for( int i = r0 + threadIdx.x; i < r1; i += BS)
        {
            const int b = A_pos[i] - base, e = A_pos[i+1] - base;
            if( beta == value_type(1))
            {
                value_type yi = y_ptr[i];
                for( int k = b; k < e; k++)
                    yi = DG_FMA( s_ax[k], s_x[k], yi);
                y_ptr[i] = yi;
            }
            else
            {
                value_type temp = 0;
                for( int k = b; k < e; k++)
                    temp = DG_FMA( s_ax[k], s_x[k], temp);
                y_ptr[i] = beta == value_type(0) ? temp : DG_FMA( beta, y_ptr[i], temp);
            }
        }
    }
    else // a single long row r0, summed in chunks of NB
    {
        value_type temp = beta == value_type(1) ? value_type(y_ptr[r0]) : value_type(0);
        for( I c = base; c < end; c += NB)
        {
            const int len = end - c < NB ? end - c : NB;
#pragma unroll
            for( int u = 0; u < NB/BS; u++)
            {
                const int k = threadIdx.x + u*BS;
                if( k < len)
                {
                    s_ax[k] = alpha*A_val[c+k];
                    s_x[k] = x_ptr[A_idx[c+k]];
                }
            }
            __syncthreads();
            if( threadIdx.x == 0)
                for( int k = 0; k < len; k++)
                    temp = DG_FMA( s_ax[k], s_x[k], temp);
            __syncthreads();
        }
        if( threadIdx.x == 0)
            y_ptr[r0] = beta == value_type(1) || beta == value_type(0) ? temp
                : DG_FMA( beta, y_ptr[r0], temp);
    }
}

//y = alpha A*x + beta y
template<class I, class V, class value_type, class C1, class C2>
void spmv_gpu_kernel(
    CSRCache_gpu& cache,
    size_t A_num_rows, size_t /*A_num_cols*/, size_t /*A_nnz*/,
    const I* A_pos , const I* A_idx, const V* A_val,
    value_type alpha, value_type beta, const C1* x_ptr, C2* y_ptr
)
{
    if( A_num_rows == 0)
        return;
    if( not cache.isUpToDate())
        cache.update( A_num_rows, A_pos);
    spmv_csr_stream_kernel<<<cache.num_blocks(), CSR_GPU_BLOCK_SIZE>>>(
        cache.row_blocks(), A_pos, A_idx, A_val, alpha, beta, x_ptr, y_ptr);
}

} // namespace detail
} // namespace dg
