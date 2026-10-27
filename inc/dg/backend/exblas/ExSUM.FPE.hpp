/*
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 *  All rights reserved.
 */

/**
 *  \file ExSUM.FPE.hpp
 *  \brief A set of routines concerning floating-point expansions
 *
 *  \authors
 *    Developers : \n
 *        Roman Iakymchuk  -- roman.iakymchuk@lip6.fr \n
 *        Sylvain Collange -- sylvain.collange@inria.fr \n
 *        Matthias Wiesenberger -- mattwi@fysik.dtu.dk
 */
#ifndef EXSUM_FPE_HPP_
#define EXSUM_FPE_HPP_
#include <algorithm> // std::fill
#include "accumulate.h"

namespace dg
{
namespace exblas
{
namespace cpu
{

/**
 * \struct FPExpansionTraits
 * \ingroup lowlevel
 * \brief This struct is meant to specify optimization or other technique used
 */
template<bool EX=false, bool FLUSHHI=false, bool H2SUM=false, bool CRF=false, bool CSWAP=false, bool B2SUM=true, bool SORT=false, bool VICT=false>
struct FPExpansionTraits
{
    static bool constexpr EarlyExit = EX;
    static bool constexpr FlushHi = FLUSHHI;
    static bool constexpr Horz2Sum = H2SUM;
    static bool constexpr CheckRangeFirst = CRF;
    static bool constexpr ConditionalSwap = CSWAP;
    static bool constexpr Biased2Sum = B2SUM;
    static bool constexpr Sort = SORT;
    static bool constexpr Victimcache = VICT;
};

/**
 * \struct FPExpansionVect
 * \ingroup lowlevel
 * \brief This struct is meant to introduce functionality for working with
 *  floating-point expansions in conjuction with superaccumulators
 */
template<typename T, int N, typename TRAITS=FPExpansionTraits<false,false> >
struct FPExpansionVect
{
    /**
     * Constructor
     * \param sa superaccumulator
     */
    FPExpansionVect(int64_t* sa);

    /**
     * This function accumulates value x to the floating-point expansion
     * \param x input value
     */
    void Accumulate(T x);

    /**
     * This function is used to flush the floating-point expansion to the superaccumulator
     */
    void Flush();

private:
    void FlushVector(T x) const;
    static T twosum(T a, T b, T & s);

    int64_t* superacc;

    // Most significant digits first!
#ifdef _MSC_VER
	_declspec(align(64)) T a[N];
#else
    T a[N] __attribute__((aligned(64)));
#endif
    T victim;
};

template<typename T, int N, typename TRAITS>
FPExpansionVect<T,N,TRAITS>::FPExpansionVect(int64_t * sa) :
    superacc(sa),
    victim(0)
{
    std::fill(a, a + N, 0);
}

template<typename T, int N, typename TRAITS> EXBLAS_UNROLL_ATTRIBUTE
void FPExpansionVect<T,N,TRAITS>::Accumulate(T x)
{
    // Experimental
    if(TRAITS::CheckRangeFirst && horizontal_or(abs(x) < abs(a[N-1]))) {
        FlushVector(x);
        return;
    }
    T s;
    for(unsigned int i = 0; i != N; ++i) {
        a[i] = twosum(a[i], x, s);
        x = s;
        if(TRAITS::EarlyExit && i != 0 && !horizontal_or(x)) return;
    }
    if(TRAITS::EarlyExit || horizontal_or(x)) {
        FlushVector(x);
    }
}

template<typename T, int N, typename TRAITS>
T FPExpansionVect<T,N,TRAITS>::twosum(T a, T b, T & s)
{
	return KnuthTwoSum(a, b, s);
}

template<typename T, int N, typename TRAITS>
void FPExpansionVect<T,N,TRAITS>::Flush()
{
    for(unsigned int i = 0; i != N; ++i)
    {
        FlushVector(a[i]);
        a[i] = 0;
    }
    if(TRAITS::Victimcache) {
        FlushVector(victim);
    }
}

template<typename T, int N, typename TRAITS> inline
void FPExpansionVect<T,N,TRAITS>::FlushVector(T x) const
{
    // TODO: update status, handle Inf/Overflow/NaN cases
    // TODO: make it work for other values of 4
    exblas::cpu::Accumulate(superacc, x);
}

}//namespace cpu
}//namespace exblas
} //namespace dg
#endif // EXSUM_FPE_HPP_
