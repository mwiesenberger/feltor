#pragma once
#include "fma.h"

namespace dg
{

// general multiply kernel
template<class real_type, class value_type>
 __global__ void ell_multiply_kernel( value_type alpha, value_type beta,
         const real_type* __restrict__  data,
         const int* __restrict__  cols_idx, const int* __restrict__  data_idx,
         const int num_rows, const int num_cols, const int blocks_per_line,
         const int n, const int size,
         const int right_size,
         const int* __restrict__  right_range,
         const value_type* __restrict__  x, value_type * __restrict__ y
         )
{
    const int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    const int grid_size = gridDim.x*blockDim.x;
    const int right_ = right_range[1]-right_range[0];
    // size counts all right_size columns but only right_ of them are computed
    const int lines = size/right_size*right_;
    //every thread takes num_rows/grid_size rows
    for( int row = thread_id; row<lines; row += grid_size)
    {
        int rr = row/right_, rrn = rr/n;
        int s=rrn/num_rows,
            i = (rrn)%num_rows,
            k = (rr)%n,
            j=right_range[0]+row%right_;
        int idx = ((s*num_rows+i)*n+k)*right_size+j;
        //idx != row ( if right_range[0] != 0)
        //y[idx]*= beta;
        // if y[I] isnan then even beta = 0 does not make it 0
        y[idx] = beta == 0 ? (value_type)0 : y[idx]*beta;
        for( int d=0; d<blocks_per_line; d++)
        {
            value_type temp=0;
            int B = (data_idx[i*blocks_per_line+d]*n+k)*n;
            int C = cols_idx[i*blocks_per_line+d];
            if( C == -1)
                continue;
            int J = (s*num_cols+C)*n;
            for( int q=0; q<n; q++) //multiplication-loop
                temp =DG_FMA( data[ B+q], x[(J+q)*right_size+j], temp);
            y[idx]=dg::detail::dg_fma( alpha, temp, y[idx]);
        }
    }

}


//specialized multiply kernel
// Every thread computes all n output rows of one cell, which means that the
// 3n input values are loaded once per cell and integer divisions by run-time
// values are done once per cell instead of once per row.
// The sum over the blocks is kept in the same order as in the general kernel
// such that both produce the same floating point results
template<class real_type, class value_type, size_t n, size_t blocks_per_line, bool unit_right>
 __global__ void ell_multiply_cell_kernel(value_type alpha, value_type beta,
         const real_type* __restrict__  data,
         const int* __restrict__  cols_idx, const int* __restrict__  data_idx,
         const int num_rows, const int num_cols,
         const int size, const int right_size,
         const int* __restrict__  right_range,
         const value_type* __restrict__  x, value_type * __restrict__ y
         )
{
    //int size = left*num_rows*right;
    const int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    const int grid_size = gridDim.x*blockDim.x;
    const int right_ = right_range[1]-right_range[0];
    const int stride = unit_right ? 1 : right_size;
    // size counts all right_size columns but only right_ of them are computed
    const int cells = size/right_size*right_;
    for( int cell = thread_id; cell<cells; cell += grid_size)
    {
        int rs = cell, j = 0;
        if( !unit_right)
        {
            rs = cell/right_;
            j = right_range[0] + cell - rs*right_;
        }
        const int s = rs/num_rows, i = rs - s*num_rows;
        value_type temp[blocks_per_line][n] = {{0}};
        for( int d=0; d<(int)blocks_per_line; d++)
        {
            const int C = cols_idx[i*blocks_per_line+d];
            if( C == -1)
                continue;
            const int B = data_idx[i*blocks_per_line+d]*n*n;
            const int J = (s*num_cols+C)*n;
            value_type xx[n];
            for( int q=0; q<(int)n; q++)
                xx[q] = x[(J+q)*stride+j];
            for( int k=0; k<(int)n; k++)
                for( int q=0; q<(int)n; q++) //multiplication-loop
                    temp[d][k] = dg::detail::dg_fma( data[ B+k*n+q], xx[q], temp[d][k]);
        }
        const int idx = (s*num_rows+i)*n*stride+j;
        for( int k=0; k<(int)n; k++)
        {
            const int I = idx + k*stride;
            value_type result = beta == 0 ? (value_type)0 : y[I]*beta;
            for( int d=0; d<(int)blocks_per_line; d++)
                result = dg::detail::dg_fma( alpha, temp[d][k], result);
            y[I] = result;
        }
    }
}

template<class real_type, class value_type, size_t n, size_t blocks_per_line>
void launch_ell_multiply_cell_kernel( size_t NUM_BLOCKS, size_t BLOCK_SIZE,
         value_type alpha, value_type beta,
         const real_type * __restrict__ data_ptr,
         const int * __restrict__ cols_ptr, const int * __restrict__ block_ptr,
         const int num_rows, const int num_cols,
         const int size, const int right_size,
         const int * __restrict__ right_range_ptr,
         const value_type * __restrict__ x_ptr, value_type * __restrict__ y_ptr)
{
    if( right_size == 1)
        ell_multiply_cell_kernel<real_type, value_type, n, blocks_per_line, true><<<NUM_BLOCKS, BLOCK_SIZE>>>
        (alpha, beta, data_ptr, cols_ptr, block_ptr, num_rows, num_cols, size,
        right_size, right_range_ptr,  x_ptr,y_ptr);
    else
        ell_multiply_cell_kernel<real_type, value_type, n, blocks_per_line, false><<<NUM_BLOCKS, BLOCK_SIZE>>>
        (alpha, beta, data_ptr, cols_ptr, block_ptr, num_rows, num_cols, size,
        right_size, right_range_ptr,  x_ptr,y_ptr);
}

template<class real_type, class value_type, size_t n>
void call_ell_multiply_kernel( value_type alpha, value_type beta,
         const real_type * __restrict__ data_ptr,
         const int * __restrict__ cols_ptr, const int * __restrict__ block_ptr,
         const int num_rows, const int num_cols, const int blocks_per_line,
         const int left_size, const int right_size,
         const int * __restrict__ right_range_ptr,
         const value_type * __restrict__ x_ptr, value_type * __restrict__ y_ptr)
{
    //set up kernel parameters
    const size_t BLOCK_SIZE = 256;
    if( blocks_per_line >= 1 && blocks_per_line <= 4)
    {
        const size_t cells = left_size*right_size*num_rows; //number of cells
        const size_t NUM_BLOCKS = std::min<size_t>((cells-1)/BLOCK_SIZE+1, 65000);
        if( blocks_per_line == 1)
            launch_ell_multiply_cell_kernel<real_type, value_type, n, 1>(NUM_BLOCKS, BLOCK_SIZE,
            alpha, beta, data_ptr, cols_ptr, block_ptr, num_rows, num_cols, cells,
            right_size, right_range_ptr,  x_ptr,y_ptr);
        else if (blocks_per_line == 2)
            launch_ell_multiply_cell_kernel<real_type, value_type, n, 2>(NUM_BLOCKS, BLOCK_SIZE,
            alpha, beta, data_ptr, cols_ptr, block_ptr, num_rows, num_cols, cells,
            right_size, right_range_ptr,  x_ptr,y_ptr);
        else if (blocks_per_line == 3)
            launch_ell_multiply_cell_kernel<real_type, value_type, n, 3>(NUM_BLOCKS, BLOCK_SIZE,
            alpha, beta, data_ptr, cols_ptr, block_ptr, num_rows, num_cols, cells,
            right_size, right_range_ptr,  x_ptr,y_ptr);
        else
            launch_ell_multiply_cell_kernel<real_type, value_type, n, 4>(NUM_BLOCKS, BLOCK_SIZE,
            alpha, beta, data_ptr, cols_ptr, block_ptr, num_rows, num_cols, cells,
            right_size, right_range_ptr,  x_ptr,y_ptr);
        return;
    }
    const size_t size = left_size*right_size*num_rows*n; //number of lines
    const size_t NUM_BLOCKS = std::min<size_t>((size-1)/BLOCK_SIZE+1, 65000);
    ell_multiply_kernel<real_type, value_type><<<NUM_BLOCKS, BLOCK_SIZE>>>
    (alpha, beta, data_ptr, cols_ptr, block_ptr, num_rows, num_cols,
    blocks_per_line, n, size, right_size, right_range_ptr,  x_ptr,y_ptr);
}


template<class real_type, template<class> class Vector>
template<class value_type>
void EllSparseBlockMat<real_type, Vector>::launch_multiply_kernel( CudaTag, value_type alpha, const value_type* x_ptr, value_type beta, value_type* y_ptr) const
{
    const real_type* data_ptr = thrust::raw_pointer_cast( &data[0]);
    const int* cols_ptr = thrust::raw_pointer_cast( &cols_idx[0]);
    const int* block_ptr = thrust::raw_pointer_cast( &data_idx[0]);
    const int* right_range_ptr = thrust::raw_pointer_cast( &right_range[0]);
    if( n == 1)
        call_ell_multiply_kernel<real_type, value_type, 1>  (alpha, beta,
            data_ptr, cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line,
            left_size, right_size, right_range_ptr,  x_ptr,y_ptr);
    else if( n == 2)
        call_ell_multiply_kernel<real_type, value_type, 2>  (alpha, beta,
            data_ptr, cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line,
            left_size, right_size, right_range_ptr,  x_ptr,y_ptr);
    else if( n == 3)
        call_ell_multiply_kernel<real_type, value_type, 3>  (alpha, beta,
            data_ptr, cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line,
            left_size, right_size, right_range_ptr,  x_ptr,y_ptr);
    else if( n == 4)
        call_ell_multiply_kernel<real_type, value_type, 4>  (alpha, beta,
            data_ptr, cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line,
            left_size, right_size, right_range_ptr,  x_ptr,y_ptr);
    else if( n == 5)
        call_ell_multiply_kernel<real_type, value_type, 5>  (alpha, beta,
            data_ptr, cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line,
            left_size, right_size, right_range_ptr,  x_ptr,y_ptr);
    else if( n == 6)
        call_ell_multiply_kernel<real_type, value_type, 6>  (alpha, beta,
            data_ptr, cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line,
            left_size, right_size, right_range_ptr,  x_ptr,y_ptr);
    else
    {
        //set up kernel parameters
        const size_t BLOCK_SIZE = 256;
        const size_t size = left_size*right_size*num_rows*n; //number of lines
        const size_t NUM_BLOCKS = std::min<size_t>((size-1)/BLOCK_SIZE+1, 65000);
        ell_multiply_kernel<real_type, value_type><<<NUM_BLOCKS, BLOCK_SIZE>>>( alpha, beta,
            data_ptr, cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line,
            n, size, right_size, right_range_ptr,  x_ptr,y_ptr);
    }
}

//////////////////// COO multiply kernel
template<class real_type, class value_type>
 __global__ void coo_multiply_kernel(
         const real_type* __restrict__  data,
         const int* __restrict__  rows_idx, const int* __restrict__  cols_idx,
         const int* __restrict__  data_idx,
         const int num_rows, const int num_cols, const int num_entries,
         const int n,
         const int left, const int right,
         value_type alpha, const value_type**  x, value_type beta,
         value_type * __restrict__ y
         )
{
    int size = left*n*right;
    const int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    const int grid_size = gridDim.x*blockDim.x;
    //every thread takes num_rows/grid_size rows
    for( int idx = thread_id; idx<size; idx += grid_size)
    {
        int s=idx/(n*right),
            k=(idx/right)%n,
            j=idx%right;
        for( int entry=0; entry<num_entries; entry++)
        {
            int I = ((s*num_rows+rows_idx[entry])*n+k)*right+j;
            value_type temp = 0;
            int B = data_idx[entry];
            int J = cols_idx[entry];
            for( int q=0; q<n; q++) //multiplication-loop
                temp = dg::detail::dg_fma( data[ (B*n + k)*n+q],
                    //x[((s*num_cols + J)*n+q)*right+j],
                    x[J][(q*left +s )*right+j],
                    temp);
            y[I] = dg::detail::dg_fma( alpha, temp, y[I]);
        }
    }
}
template<class real_type, class value_type, int n>
 __global__ void coo_multiply_kernel(
         const real_type* __restrict__  data,
         const int* __restrict__  rows_idx, const int* __restrict__  cols_idx,
         const int* __restrict__  data_idx,
         const int num_rows, const int num_cols, const int num_entries,
         const int left, const int right,
         value_type alpha, const value_type**  x, value_type beta,
         value_type * __restrict__ y
         )
{
    int size = left*n*right;
    const int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    const int grid_size = gridDim.x*blockDim.x;
    //every thread takes num_rows/grid_size rows
    for( int idx = thread_id; idx<size; idx += grid_size)
    {
        int s=idx/(n*right),
            k=(idx/right)%n,
            j=idx%right;
        for( int entry=0; entry<num_entries; entry++)
        {
            int I = ((s*num_rows+rows_idx[entry])*n+k)*right+j;
            value_type temp = 0;
            int B = data_idx[entry];
            int J = cols_idx[entry];
            for( int q=0; q<n; q++) //multiplication-loop
                temp = dg::detail::dg_fma( data[ (B*n + k)*n+q],
                    //x[((s*num_cols + J)*n+q)*right+j],
                    x[J][(q*left +s )*right+j],
                    temp);
            y[I] = dg::detail::dg_fma( alpha, temp, y[I]);
        }
    }
}

template<class real_type, template<class> class Vector>
template<class value_type>
void CooSparseBlockMat<real_type, Vector>::launch_multiply_kernel( CudaTag, value_type alpha, const value_type** x_ptr, value_type beta, value_type* y_ptr) const
{
    if( num_entries == 0)
        return;
    assert( beta == 1 && "Beta != 1 yields wrong results in CooSparseBlockMat!!");
    //set up kernel parameters
    const size_t BLOCK_SIZE = 256;
    const size_t size = left_size*right_size*n;
    const size_t NUM_BLOCKS = std::min<size_t>((size-1)/BLOCK_SIZE+1, 65000);

    const real_type* data_ptr = thrust::raw_pointer_cast( data.data());
    const int* rows_ptr = thrust::raw_pointer_cast( rows_idx.data());
    const int* cols_ptr = thrust::raw_pointer_cast( cols_idx.data());
    const int* block_ptr = thrust::raw_pointer_cast( data_idx.data());
    if( n == 1)
        coo_multiply_kernel<real_type, value_type, 1> <<<NUM_BLOCKS, BLOCK_SIZE>>> (
            data_ptr, rows_ptr, cols_ptr, block_ptr, num_rows, num_cols,
            num_entries, left_size, right_size, alpha, x_ptr, beta, y_ptr);
    else if( n == 2)
        coo_multiply_kernel<real_type, value_type, 2> <<<NUM_BLOCKS, BLOCK_SIZE>>> (
            data_ptr, rows_ptr, cols_ptr, block_ptr, num_rows, num_cols,
            num_entries, left_size, right_size, alpha, x_ptr, beta, y_ptr);
    else if( n == 3)
        coo_multiply_kernel<real_type, value_type, 3> <<<NUM_BLOCKS, BLOCK_SIZE>>> (
            data_ptr, rows_ptr, cols_ptr, block_ptr, num_rows, num_cols,
            num_entries, left_size, right_size, alpha, x_ptr, beta, y_ptr);
    else if( n == 4)
        coo_multiply_kernel<real_type, value_type, 4> <<<NUM_BLOCKS, BLOCK_SIZE>>> (
            data_ptr, rows_ptr, cols_ptr, block_ptr, num_rows, num_cols,
            num_entries, left_size, right_size, alpha, x_ptr, beta, y_ptr);
    else
        coo_multiply_kernel<real_type, value_type> <<<NUM_BLOCKS, BLOCK_SIZE>>> (
            data_ptr, rows_ptr, cols_ptr, block_ptr, num_rows, num_cols,
            num_entries, n, left_size, right_size, alpha, x_ptr, beta, y_ptr);
}

}//namespace dg
