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
#include "exdot_serial.h"
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
    AddSuperaccs( acc1, acc2, 1);
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
    using value_type = get_result_type<Functor, PointerOrValues...>;
    constexpr unsigned K = ValueTraits<value_type>::K; // one superaccumulator per component
    int const linesize = 16;    // * sizeof(int32_t)
    int maxthreads = omp_get_max_threads();
    // per component: the superaccumulators of all threads and their progress
    std::array<std::vector<int64_t>, K> acc;
    std::array<std::vector<std::atomic<int32_t>>, K> ready;
    for( unsigned k=0; k<K; k++)
    {
        acc[k].assign( maxthreads*BIN_COUNT, 0);
        ready[k] = std::vector<std::atomic<int32_t>>( maxthreads * linesize); // value-initialized to 0
    }
    std::vector<int> error( maxthreads, false); // not vector<bool>: threads write concurrently

    // Each thread calls its own copy of f
    #pragma omp parallel firstprivate(f)
    {
        unsigned int tid = omp_get_thread_num();
        unsigned int tnum = omp_get_num_threads();
        auto thread_acc = [&acc, tid]( size_t k){ return &acc[k][tid*BIN_COUNT];};

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
            for( int j=0; j<8; j++)
            {
                value_type v = f( get_element(xs,i+j)...);
                for( unsigned k=0; k<K; k++)
                    tmp[k][j] = get_result( v, k);
            }
            for( unsigned k=0; k<K; k++)
            {
                vcl::Vec8d x  = vcl::Vec8d().load( tmp[k]);
                vcl::Vec8db finite = vcl::is_finite( x);
                if( !vcl::horizontal_and( finite) ) error[tid] = true;
                cache[k].Accumulate(x);
            }
        }
        if( tid+1==tnum && r != N-1) {
            r+=1;
            //accumulate remainder
            double tmp[K][8];
            for( int j=0; j<N-r; j++)
            {
                value_type v = f( get_element(xs,r+j)...);
                for( unsigned k=0; k<K; k++)
                    tmp[k][j] = get_result( v, k);
            }
            for( unsigned k=0; k<K; k++)
            {
                vcl::Vec8d x  = vcl::Vec8d().load_partial( N-r, tmp[k]);
                vcl::Vec8db finite = vcl::is_finite( x);
                if( !vcl::horizontal_and( finite) ) error[tid] = true;
                cache[k].Accumulate(x);
            }
        }
#else// EXBLAS_WITHOUT_VCL
        using Cache = FPExpansionVect<double, NBFPE, FPExpansionTraits<true> >;
        auto cache = make_caches<Cache>( thread_acc, std::make_index_sequence<K>());
        int l = ((tid * int64_t(N)) / tnum);
        int r = ((((tid+1) * int64_t(N)) / tnum) ) - 1;
        for(int i = l; i <= r; i++) {
            value_type v = f( get_element(xs,i)...);
            for( unsigned k=0; k<K; k++)
            {
                double x = get_result( v, k);
                if( !std::isfinite(x) ) error[tid] = true;
                cache[k].Accumulate(x);
            }
        }
#endif// EXBLAS_WITHOUT_VCL
        for( unsigned k=0; k<K; k++)
        {
            cache[k].Flush();
            int imin=IMIN, imax=IMAX;
            Normalize(thread_acc(k), imin, imax);
        }
        for( unsigned k=0; k<K; k++)
            Reduction(tid, tnum, ready[k], acc[k], linesize);
    }
    for( unsigned k=0; k<K; k++)
    for( int i=IMIN; i<=IMAX; i++)
        h_superacc[k*BIN_COUNT+i] = acc[k][i];
    for ( int i=0; i<maxthreads; i++)
        if( error[i] == true) *err = true;
}
}//namespace cpu
///@endcond

/*! @brief The size below which \c exsum_omp calls the serial \c exsum_cpu
 *
 * An exact sum is cheap compared to the parallel region and the tree
 * reduction of the superaccumulators. Serial is faster up to at least this
 * size (measured on a 6 core i5-14400F only; the crossover is at about
 * 8000-16000 there and depends on the machine)
 */
constexpr unsigned MIN_OMP_SIZE = 4096;

///@brief OpenMP parallel version of exact sum
///@copydoc hide_exsum
///@copydoc hide_hostacc
///@note Sizes smaller than \c exblas::MIN_OMP_SIZE are summed serially with \c exsum_cpu
template<size_t NBFPE=8, class Functor, class ...PointerOrValues>
void exsum_omp(unsigned size, int64_t* h_superacc, int* status, Functor f, PointerOrValues ...xs){
    static_assert( ValueTraits<get_result_type<Functor, PointerOrValues...>>::K > 0,
        "The values of f must be float, double or complex numbers of these");
    if( size < MIN_OMP_SIZE)
    {
        exsum_cpu<NBFPE>( size, h_superacc, status, f, xs...);
        return;
    }
    bool error = false;
    cpu::ExSUMFPE<NBFPE>( (int)size, h_superacc, &error, f, xs...);
    *status = 0;
    if( error ) *status = 1;
}

}//namespace exblas
} //namespace dg
