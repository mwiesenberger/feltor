#include <omp.h>
#include "config.h"
#include "fma.h"

//for the fmas it is important to activate -mfma compiler flag

namespace dg{

// general multiply kernel
template<class real_type, class value_type>
void ell_omp_multiply_kernel( value_type alpha, value_type beta,
         const real_type * RESTRICT data, const int * RESTRICT cols_idx,
         const int * RESTRICT data_idx,
         const int num_rows, const int num_cols, const int blocks_per_line,
         const int n,
         const int left_size, const int right_size,
         const int * RESTRICT right_range,
         const value_type * RESTRICT x, value_type * RESTRICT y
         )
{
    std::vector<int> J(blocks_per_line), B(blocks_per_line);
#pragma omp for nowait //manual collapse(2)
	for( int si = 0; si<left_size*num_rows; si++)
	{
		int s = si / num_rows;
		int i = si % num_rows;
        for( int d=0; d<blocks_per_line; d++)
        {
            int C = cols_idx[i*blocks_per_line+d];
            J[d] = ( C  == -1 ? -1 : (s*num_cols+C)*n);
        }
        for( int k=0; k<n; k++)
        {
            for( int d=0; d<blocks_per_line; d++)
                B[d] = (data_idx[i*blocks_per_line+d]*n+k)*n;
            for( int j=right_range[0]; j<right_range[1]; j++)
            {
                int I = ((s*num_rows + i)*n+k)*right_size+j;
                // if y[I] isnan then even beta = 0 does not make it 0
                y[I] = beta == 0 ? (value_type)0 : y[I]*beta;
                for( int d=0; d<blocks_per_line; d++)
                {
                    value_type temp = 0;
                    if( J[d] == -1)
                        continue;
                    for( int q=0; q<n; q++) //multiplication-loop
                        temp = DG_FMA(data[ B[d]+q],
                                x[(J[d]+q)*right_size+j],
                                temp);
                    y[I] = DG_FMA(alpha, temp, y[I]);
                }
            }
        }
    }
}
//specialized multiply kernel
template<class real_type, class value_type, int n, int blocks_per_line>
void ell_omp_multiply_kernel( value_type alpha, value_type beta,
         const real_type * RESTRICT data, const int * RESTRICT cols_idx,
         const int * RESTRICT data_idx,
         const int num_rows, const int num_cols,
         const int left_size, const int right_size,
         const int * RESTRICT right_range,
         const value_type * RESTRICT x, value_type * RESTRICT y
         )
{
    //basically we check which direction is the largest and parallelize that one
    if(right_size==1)
    {
    // trivial means that the data blocks do not change among rows
    // that are not the first or the last row (where BCs usually live)
    bool trivial = true;
    if( num_rows < 4) // need at least 3 rows for this to make sense
        trivial = false;
    for( int i=2; i<num_rows-1; i++)
        for( int d=0; d<blocks_per_line; d++)
        {
            if( data_idx[i*blocks_per_line+d]
                    != data_idx[blocks_per_line+d]) trivial = false;
        }
    if(trivial)
    {
    value_type xprivate[blocks_per_line*n];
    real_type dprivate[blocks_per_line*n*n];
    for( int d=0; d<blocks_per_line; d++)
    for( int k=0; k<n; k++)
    for( int q=0; q<n; q++)
    {
        int B = data_idx[blocks_per_line+d];
        dprivate[(k*blocks_per_line+d)*n+q] = data[(B*n+k)*n+q];
    }
    #pragma omp for nowait
    for( int s=0; s<left_size; s++)
    {
        for( int i=0; i<1; i++)
        {
            for( int d=0; d<blocks_per_line; d++)
            {
                int C = cols_idx[i*blocks_per_line+d];
                int J = (s*num_cols+C)*n;
                for(int q=0; q<n; q++)
                    xprivate[d*n+q] = (C == -1 ? 0 : x[J+q]);
            }
            for( int k=0; k<n; k++)
            {
                value_type temp[blocks_per_line] = {0};
                for( int d=0; d<blocks_per_line; d++)
                {
                    int B = (data_idx[i*blocks_per_line+d]*n+k)*n;
                    for( int q=0; q<n; q++) //multiplication-loop
                        temp[d] = DG_FMA(data[B+q], xprivate[d*n+q], temp[d]);
                }
                int I = ((s*num_rows + i)*n+k);
                // if y[I] isnan then even beta = 0 does not make it 0
                y[I] = beta == 0 ? (value_type)0 : y[I]*beta;
                for( int d=0; d<blocks_per_line; d++)
                    y[I] = DG_FMA(alpha, temp[d], y[I]);
            }
        }
        #if _OPENMP > 201300 // OpenMP is >= v4.0
        #pragma omp SIMD //very important for KNL (6.8.25 : Clang does not vectorize this)
        #endif
        for( int i=1; i<num_rows-1; i++)
        {
            for( int k=0; k<n; k++)
            {
                int I = ((s*num_rows + i)*n+k);
                // if y[I] isnan then even beta = 0 does not make it 0
                y[I] = beta == 0 ? (value_type)0 : y[I]*beta;
                int B = n*blocks_per_line*k;
                for( int d=0; d<blocks_per_line; d++)
                {
                    value_type temp = 0;
                    int C = cols_idx[i*blocks_per_line+d];
                    if( C == -1)
                        continue;
                    for( int q=0; q<n; q++)
                    {
                        int J = (s*num_cols+C)*n+q;
                        temp = DG_FMA( dprivate[B+d*n+q], x[J], temp);
                    }
                    y[I] = DG_FMA(alpha, temp, y[I]);
                }
            }
        }
        for( int i=num_rows-1; i<num_rows; i++)
        {
            for( int d=0; d<blocks_per_line; d++)
            {
                int C = cols_idx[i*blocks_per_line+d];
                int J = (s*num_cols+C)*n;
                for(int q=0; q<n; q++)
                    xprivate[d*n+q] = (C == -1 ? 0 : x[J+q]);
            }
            for( int k=0; k<n; k++)
            {
                value_type temp[blocks_per_line] = {0};
                for( int d=0; d<blocks_per_line; d++)
                {
                    int B = (data_idx[i*blocks_per_line+d]*n+k)*n;
                    for( int q=0; q<n; q++) //multiplication-loop
                        temp[d] = DG_FMA( data[B+q], xprivate[d*n+q], temp[d]);
                }
                int I = ((s*num_rows + i)*n+k);
                // if y[I] isnan then even beta = 0 does not make it 0
                y[I] = beta == 0 ? (value_type)0 : y[I]*beta;
                for( int d=0; d<blocks_per_line; d++)
                    y[I] = DG_FMA(alpha, temp[d], y[I]);
            }
        }
    }
    } //trivial
    else // not trivial
    {
    value_type xprivate[blocks_per_line*n];
    #pragma omp for nowait
    for( int s=0; s<left_size; s++)
    for( int i=0; i<num_rows; i++)
    {
        for( int d=0; d<blocks_per_line; d++)
        {
            int C = cols_idx[i*blocks_per_line+d];
            int J = (s*num_cols+C)*n;
            for(int q=0; q<n; q++)
                xprivate[d*n+q] = (C == -1 ? 0 : x[J+q]);
        }
        for( int k=0; k<n; k++)
        {
            value_type temp[blocks_per_line] = {0};
            for( int d=0; d<blocks_per_line; d++)
            {
                int B = (data_idx[i*blocks_per_line+d]*n+k)*n;
                for( int q=0; q<n; q++) //multiplication-loop
                    temp[d] = DG_FMA( data[B+q], xprivate[d*n+q], temp[d]);
            }
            int I = ((s*num_rows + i)*n+k);
            // if y[I] isnan then even beta = 0 does not make it 0
            y[I] = beta == 0 ? (value_type)0 : y[I]*beta;
            for( int d=0; d<blocks_per_line; d++)
                y[I] = DG_FMA(alpha, temp[d], y[I]);
        }
    }
    }// trivial
    }// right_size==1
    else // right_size != 1
    {
    real_type dprivate[blocks_per_line*n];
    // J and Jd are long instead of int on purpose: with int, gcc 13 does not
    // vectorize the j loops below for some template instantiations, and
    // derivatives in y are up to 1.4x slower. Widening only I or j does not help.
    long J[blocks_per_line];
    if( !( (right_range[1]-right_range[0]) > 100*left_size*num_rows*n )) //typically a derivative in y ( Ny*Nz >~ Nx)
    {
        #pragma omp for nowait
        for (int sik = 0; sik < left_size*num_rows*n; sik++)
        {
            int s = sik / (num_rows*n);
            int i = (sik % (num_rows*n)) / n;
            int k = (sik % (num_rows*n)) % n;

            for( int d=0; d<blocks_per_line; d++)
            {
                int C = cols_idx[i*blocks_per_line+d];
                J[d] = ( C == -1 ? -1 :(s*num_cols+C)*n );
                int B = (data_idx[i*blocks_per_line+d]*n+k)*n;
                for(int q=0; q<n; q++)
                    dprivate[d*n+q] = data[B+q];
            }
            #if _OPENMP > 201300 // OpenMP is >= v4.0
            #pragma omp SIMD //very important for KNL
            #endif
            for( int j=right_range[0]; j<right_range[1]; j++)
            {
                int I = ((s*num_rows + i)*n+k)*right_size+j;
                // if y[I] isnan then even beta = 0 does not make it 0
                y[I] = beta == 0 ? (value_type)0 : y[I]*beta;
                for( int d=0; d<blocks_per_line; d++)
                {
                    long Jd = J[d];
                    if ( Jd == -1)
                        continue;
                    value_type temp = 0;
                    for( int q=0; q<n; q++) //multiplication-loop
                        temp = DG_FMA( dprivate[ d*n+q],
                                    x[(Jd+q)*right_size+j],
                                    temp);
                    y[I] = DG_FMA(alpha, temp, y[I]);
                }
            }
        }
    }
    else //typically a derivative in z (since n*n*Nx*Ny > 100*Nz)
    {

        for (int sik = 0; sik < left_size*num_rows*n; sik++)
        {
            int s = sik / (num_rows*n);
            int i = (sik % (num_rows*n)) / n;
            int k = (sik % (num_rows*n)) % n;

            for( int d=0; d<blocks_per_line; d++)
            {
                int C = cols_idx[i*blocks_per_line+d];
                J[d] = ( C == -1 ? -1 :(s*num_cols+C)*n );
                int B = (data_idx[i*blocks_per_line+d]*n+k)*n;
                for(int q=0; q<n; q++)
                    dprivate[d*n+q] = data[B+q];
            }
            #pragma omp for SIMD nowait
            for( int j=right_range[0]; j<right_range[1]; j++)
            {
                int I = ((s*num_rows + i)*n+k)*right_size+j;
                // if y[I] isnan then even beta = 0 does not make it 0
                y[I] = beta == 0 ? (value_type)0 : y[I]*beta;
                for( int d=0; d<blocks_per_line; d++)
                {
                    long Jd = J[d];
                    if( Jd == -1)
                        continue;
                    value_type temp = 0;
                    for( int q=0; q<n; q++) //multiplication-loop
                        temp = DG_FMA( dprivate[ d*n+q],
                                    x[(Jd+q)*right_size+j],
                                    temp);
                    y[I] = DG_FMA(alpha, temp, y[I]);
                }
            }
        }
        }
    }
}

template<class real_type, class value_type, int n>
void call_ell_omp_multiply_kernel( value_type alpha, value_type beta,
         const real_type * RESTRICT data_ptr, const int * RESTRICT cols_ptr,
         const int * RESTRICT block_ptr,
         const int num_rows, const int num_cols, const int blocks_per_line,
         const int left_size, const int right_size,
         const int * RESTRICT right_range_ptr,
         const value_type * RESTRICT x_ptr, value_type * RESTRICT y_ptr)
{
    if( blocks_per_line == 1)
        ell_omp_multiply_kernel<real_type, value_type, n, 1>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, left_size, right_size,
        right_range_ptr,  x_ptr,y_ptr);
    else if (blocks_per_line == 2)
        ell_omp_multiply_kernel<real_type, value_type, n, 2>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, left_size, right_size,
        right_range_ptr,  x_ptr,y_ptr);
    else if (blocks_per_line == 3)
        ell_omp_multiply_kernel<real_type, value_type, n, 3>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, left_size, right_size,
        right_range_ptr,  x_ptr,y_ptr);
    else if (blocks_per_line == 4)
        ell_omp_multiply_kernel<real_type, value_type, n, 4>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, left_size, right_size,
        right_range_ptr,  x_ptr,y_ptr);
    else
        ell_omp_multiply_kernel<real_type, value_type>  (alpha, beta, data_ptr, cols_ptr,
        block_ptr, num_rows, num_cols, blocks_per_line, n, left_size,
        right_size, right_range_ptr,  x_ptr,y_ptr);
}


template<class real_type, template<class> class Vector>
template<class value_type>
void EllSparseBlockMat<real_type, Vector>::launch_multiply_kernel( OmpTag, value_type alpha, const value_type* x_ptr, value_type beta, value_type* y_ptr) const
{
    const real_type* data_ptr = thrust::raw_pointer_cast( &data[0]);
    const int* cols_ptr = thrust::raw_pointer_cast( &cols_idx[0]);
    const int* block_ptr = thrust::raw_pointer_cast( &data_idx[0]);
    const int* right_range_ptr = thrust::raw_pointer_cast( &right_range[0]);
    if( n == 1)
        call_ell_omp_multiply_kernel<real_type, value_type, 1>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line, left_size,
        right_size, right_range_ptr,  x_ptr,y_ptr);

    else if( n == 2)
        call_ell_omp_multiply_kernel<real_type, value_type, 2>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line, left_size,
        right_size, right_range_ptr,  x_ptr,y_ptr);
    else if( n == 3)
        call_ell_omp_multiply_kernel<real_type, value_type, 3>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line, left_size,
        right_size, right_range_ptr,  x_ptr,y_ptr);
    else if( n == 4)
        call_ell_omp_multiply_kernel<real_type, value_type, 4>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line, left_size,
        right_size, right_range_ptr,  x_ptr,y_ptr);
    else if( n == 5)
        call_ell_omp_multiply_kernel<real_type, value_type, 5>  (alpha, beta, data_ptr,
        cols_ptr, block_ptr, num_rows, num_cols, blocks_per_line, left_size,
        right_size, right_range_ptr,  x_ptr,y_ptr);
    else
        ell_omp_multiply_kernel<real_type, value_type> ( alpha, beta, data_ptr, cols_ptr,
        block_ptr, num_rows, num_cols, blocks_per_line, n, left_size,
        right_size, right_range_ptr,  x_ptr,y_ptr);
}

template<class real_type, class value_type, template<class> class Vector>
void coo_omp_multiply_kernel( value_type alpha, const value_type** x, value_type /*beta*/, value_type* RESTRICT y, const CooSparseBlockMat<real_type, Vector>& m )
{
    // One loop per entry is faster than all entries in one loop (hypothesis:
    // the loop over sj vectorises). A static schedule without simd gives
    // every thread the same sj in every loop, so there is no race on y and
    // each y[I] sums its entries in the same order as before.
    // Entry values are read before the loop over sj: gcc cannot tell that the
    // writes to y leave them unchanged, and with them inside the loop gcc 13
    // warns about undefined behaviour (-Waggressive-loop-optimizations) for n=4
    const int n = m.n, left_size = m.left_size, right_size = m.right_size;
    const int num_rows = m.num_rows, size = left_size*right_size;
    const real_type* RESTRICT data = thrust::raw_pointer_cast( m.data.data());
    for (int i = 0; i < m.num_entries; i++)
    {
        const int row = m.rows_idx[i], B = m.data_idx[i];
        const value_type* RESTRICT xi = x[m.cols_idx[i]];
        #pragma omp for schedule(static) nowait
        for (int sj = 0; sj < size; sj++)
        {
            int s = sj / right_size;
            int j = sj % right_size;
            for( int k=0; k<n; k++)
            {
                int I = ((s*num_rows + row)*n + k)*right_size + j;
                value_type temp = 0;
                for (int q = 0; q < n; q++) //multiplication-loop
                    temp = DG_FMA(data[(B * n + k)*n + q],
                        //x[((s*m.num_cols + m.cols_idx[i])*n+q)*m.right_size+j],
                        xi[q*size +sj],
                        temp);
                y[I] = DG_FMA(alpha, temp, y[I]);
            }
        }
    }
}
template<class real_type, class value_type, int n, template<class > class Vector>
void coo_omp_multiply_kernel( value_type alpha, const value_type** x, value_type /*beta*/, value_type* RESTRICT y, const CooSparseBlockMat<real_type, Vector>& m )
{
    // see the general kernel for the loop order, the schedule and the locals
    const int left_size = m.left_size, right_size = m.right_size;
    const int num_rows = m.num_rows, size = left_size*right_size;
    const real_type* RESTRICT data = thrust::raw_pointer_cast( m.data.data());
    for (int i = 0; i < m.num_entries; i++)
    {
        const int row = m.rows_idx[i], B = m.data_idx[i];
        const value_type* RESTRICT xi = x[m.cols_idx[i]];
        #pragma omp for schedule(static) nowait
        for (int sj = 0; sj < size; sj++)
        {
            int s = sj / right_size;
            int j = sj % right_size;
            for( int k=0; k<n; k++)
            {
                int I = ((s*num_rows + row)*n + k)*right_size + j;
                value_type temp = 0;
                for (int q = 0; q < n; q++) //multiplication-loop
                    temp = DG_FMA(data[(B * n + k)*n + q],
                        //x[((s*m.num_cols + m.cols_idx[i])*n+q)*m.right_size+j],
                        xi[q*size +sj],
                        temp);
                y[I] = DG_FMA(alpha, temp, y[I]);
            }
        }
    }
}
template<class real_type, template<class> class Vector>
template<class value_type>
void CooSparseBlockMat<real_type, Vector>::launch_multiply_kernel( OmpTag, value_type alpha, const value_type** x, value_type beta, value_type* RESTRICT y) const
{
    if( num_entries==0)
        return;
    assert( beta == 1 && "Beta != 1 yields wrong results in CooSparseBlockMat!!");
    if( n == 1)
        coo_omp_multiply_kernel<real_type, value_type, 1>( alpha, x, beta, y, *this);
    else if( n == 2)
        coo_omp_multiply_kernel<real_type, value_type, 2>( alpha, x, beta, y, *this);
    else if( n == 3)
        coo_omp_multiply_kernel<real_type, value_type, 3>( alpha, x, beta, y, *this);
    else if( n == 4)
        coo_omp_multiply_kernel<real_type, value_type, 4>( alpha, x, beta, y, *this);
    else if( n == 5)
        coo_omp_multiply_kernel<real_type, value_type, 5>( alpha, x, beta, y, *this);
    else
        coo_omp_multiply_kernel<real_type, value_type>( alpha, x, beta, y, *this);
}

}//namespace dg
