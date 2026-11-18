/*
 * %%%%%%%%%%%%%%%%%%%%%%%Original development%%%%%%%%%%%%%%%%%%%%%%%%%
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 * %%%%%%%%%%%%%%%%%%%%%%%Modifications and further additions%%%%%%%%%%
 *  Matthias Wiesenberger, 2017, within FELTOR and EXBLAS licenses
 */

/**
 *  @file exdot_omp.h
 *  @brief OpenMP version of exdot
 *
 *  @authors
 *    Developers : \n
 *        Roman Iakymchuk  -- roman.iakymchuk@lip6.fr \n
 *        Sylvain Collange -- sylvain.collange@inria.fr \n
 *        Matthias Wiesenberger -- mattwi@fysik.dtu.dk
 */
#pragma once
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <atomic>

#include "accumulate.h"
#include "ExSUM.FPE.hpp"
#include <omp.h>

namespace dg
{
namespace exblas{
///@cond
namespace cpu{

/**
 * \brief Parallel reduction step
 *
 * \param step step among threads
 * \param acc1 superaccumulator of the first thread
 * \param acc2 superaccumulator of the second thread
 * \param ready level the second thread has reached
 */
inline void ReductionStep(int step, int64_t * acc1, int64_t * acc2,
    const std::atomic<int32_t> * ready)
{
    // Wait for thread 2 to be ready; acquire makes its accumulator visible
    while(ready->load( std::memory_order_acquire) < step) {
#ifndef EXBLAS_WITHOUT_VCL
        _mm_pause();
#endif//EXBLAS_WITHOUT_VCL
    }
    int imin = IMIN, imax = IMAX;
    Normalize( acc1, imin, imax);
    imin = IMIN, imax = IMAX;
    Normalize( acc2, imin, imax);
    for(int i = IMIN; i <= IMAX; ++i) {
        acc1[i] += acc2[i];
    }
}

/**
 * \brief Final step of summation -- Parallel reduction among threads
 *
 * \param tid thread ID
 * \param tnum number of threads
 * \param acc superaccumulator
 */
inline void Reduction(unsigned int tid, unsigned int tnum, std::vector<std::atomic<int32_t>>& ready,
    std::vector<int64_t>& acc, int const linesize)
{
    // Custom tree reduction
    for(unsigned int s = 1; (unsigned)(1 << (s-1)) < tnum; ++s)
    {
        // 1<<(s-1) = 0001, 0010, 0100, ... = 1,2,4,8,16,...
        // set: ready for level s; release publishes our accumulator
        // (a volatile flag is neither atomic nor a memory fence)
        ready[tid * linesize].store( s, std::memory_order_release);
        if(tid % (1 << s) == 0) { //1<<s = 2,4,8,16,32,...
            //only the tid thread executes this block, tid2 just sets ready
            unsigned int tid2 = tid | (1 << (s-1)); //effectively adds 1, 2, 4,...
            if(tid2 < tnum) {
                ReductionStep(s, &acc[tid*BIN_COUNT], &acc[tid2*BIN_COUNT],
                    &ready[tid2 * linesize]);
            }
        }
    }
}

// Exact sum of f( x_0[i], x_1[i], ...) into h_superacc
// Each thread vectorises like ExSUMFPE_cpu (see there why that matters)
template<size_t NBFPE, class Functor, class ...PointerOrValues>
void ExSUMFPE(int N, int64_t* h_superacc, bool* err, Functor f, PointerOrValues ...xs) {
    using value_type = return_type<Functor, PointerOrValues...>;
    constexpr unsigned K = ValueTraits<value_type>::K; // one superaccumulator per component
    int const linesize = 16;    // * sizeof(int32_t)
    int maxthreads = omp_get_max_threads();
    // per component: the superaccumulators of all threads and their progress
    std::array<std::vector<int64_t>, K> acc;
    std::array<std::vector<std::atomic<int32_t>>, K> ready;
    for( unsigned c=0; c<K; c++)
    {
        acc[c].assign( maxthreads*BIN_COUNT, 0);
        ready[c] = std::vector<std::atomic<int32_t>>( maxthreads * linesize); // value-initialized to 0
    }
    std::vector<int> error( maxthreads, false); // not vector<bool>: threads write concurrently

    // Each thread calls its own copy of f
    #pragma omp parallel firstprivate(f)
    {
        unsigned int tid = omp_get_thread_num();
        unsigned int tnum = omp_get_num_threads();
        auto thread_acc = [&acc, tid]( size_t c){ return &acc[c][tid*BIN_COUNT];};

#ifndef EXBLAS_WITHOUT_VCL
        using Cache = FPExpansionVect<vcl::Vec8d, NBFPE, FPExpansionTraits<true> >;
        auto cache = make_caches<Cache>( thread_acc, std::make_index_sequence<K>());
        int l = ((tid * int64_t(N)) / tnum) & ~7ul; // & ~7ul == round down to multiple of 8
        int r = ((((tid+1) * int64_t(N)) / tnum) & ~7ul) - 1;

        for(int i = l; i < r; i+=8) {
#ifndef _MSC_VER
            asm ("# myloop");
#endif
            double tmp[K][8];
            for( int k=0; k<8; k++)
            {
                value_type v = f( get_element(xs,i+k)...);
                for( unsigned c=0; c<K; c++)
                    tmp[c][k] = ValueTraits<value_type>::get( v, c);
            }
            for( unsigned c=0; c<K; c++)
            {
                vcl::Vec8d x  = vcl::Vec8d().load( tmp[c]);
                vcl::Vec8db finite = vcl::is_finite( x);
                if( !vcl::horizontal_and( finite) ) error[tid] = true;
                cache[c].Accumulate(x);
            }
        }
        if( tid+1==tnum && r != N-1) {
            r+=1;
            //accumulate remainder
            double tmp[K][8];
            for( int k=0; k<N-r; k++)
            {
                value_type v = f( get_element(xs,r+k)...);
                for( unsigned c=0; c<K; c++)
                    tmp[c][k] = ValueTraits<value_type>::get( v, c);
            }
            for( unsigned c=0; c<K; c++)
            {
                vcl::Vec8d x  = vcl::Vec8d().load_partial( N-r, tmp[c]);
                vcl::Vec8db finite = vcl::is_finite( x);
                if( !vcl::horizontal_and( finite) ) error[tid] = true;
                cache[c].Accumulate(x);
            }
        }
#else// EXBLAS_WITHOUT_VCL
        using Cache = FPExpansionVect<double, NBFPE, FPExpansionTraits<true> >;
        auto cache = make_caches<Cache>( thread_acc, std::make_index_sequence<K>());
        int l = ((tid * int64_t(N)) / tnum);
        int r = ((((tid+1) * int64_t(N)) / tnum) ) - 1;
        for(int i = l; i <= r; i++) {
            value_type v = f( get_element(xs,i)...);
            for( unsigned c=0; c<K; c++)
            {
                double x = ValueTraits<value_type>::get( v, c);
                if( !std::isfinite(x) ) error[tid] = true;
                cache[c].Accumulate(x);
            }
        }
#endif// EXBLAS_WITHOUT_VCL
        for( unsigned c=0; c<K; c++)
        {
            cache[c].Flush();
            int imin=IMIN, imax=IMAX;
            Normalize(thread_acc(c), imin, imax);
        }
        for( unsigned c=0; c<K; c++)
            Reduction(tid, tnum, ready[c], acc[c], linesize);
    }
    for( unsigned c=0; c<K; c++)
    for( int i=IMIN; i<=IMAX; i++)
        h_superacc[c*BIN_COUNT+i] = acc[c][i];
    for ( int i=0; i<maxthreads; i++)
        if( error[i] == true) *err = true;
}
}//namespace cpu
///@endcond

///@brief OpenMP parallel version of exact sum
///@copydoc hide_exsum
///@copydoc hide_hostacc
template<size_t NBFPE=8, class Functor, class ...PointerOrValues>
void exsum_omp(unsigned size, int64_t* h_superacc, int* status, Functor f, PointerOrValues ...xs){
    bool error = false;
    cpu::ExSUMFPE<NBFPE>( (int)size, h_superacc, &error, f, xs...);
    *status = 0;
    if( error ) *status = 1;
}

}//namespace exblas
} //namespace dg
