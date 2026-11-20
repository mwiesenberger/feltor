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
#include <array>
#include <utility> // std::index_sequence
#include "accumulate.h"

namespace dg
{
namespace exblas
{
namespace cpu
{

/**
 * \struct FPExpansionVect
 * \ingroup lowlevel
 * \brief This struct is meant to introduce functionality for working with
 *  floating-point expansions in conjuction with superaccumulators
 */
template<typename T, int N>
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

    int64_t* superacc;

    // Most significant digits first!
#ifdef _MSC_VER
	_declspec(align(64)) T a[N];
#else
    T a[N] __attribute__((aligned(64)));
#endif
};

template<typename T, int N>
FPExpansionVect<T,N>::FPExpansionVect(int64_t * sa) :
    superacc(sa)
{
    std::fill(a, a + N, 0);
}

template<typename T, int N>
void FPExpansionVect<T,N>::Accumulate(T x)
{
    T s;
    for(unsigned int i = 0; i != N; ++i) {
        a[i] = KnuthTwoSum(a[i], x, s);
        x = s;
        // Early exit once the rounding error is zero
        if(i != 0 && !horizontal_or(x)) return;
    }
    FlushVector(x);
}

template<typename T, int N>
void FPExpansionVect<T,N>::Flush()
{
    for(unsigned int i = 0; i != N; ++i)
    {
        FlushVector(a[i]);
        a[i] = 0;
    }
}

template<typename T, int N> inline
void FPExpansionVect<T,N>::FlushVector(T x) const
{
    exblas::cpu::Accumulate(superacc, x);
}

// One cache per component of a value, the I-th one flushes to acc(I)
template<class Cache, class AccPointer, size_t ...I>
std::array<Cache, sizeof...(I)> make_caches( AccPointer acc, std::index_sequence<I...>)
{
    return {{ Cache( acc(I))...}};
}

}//namespace cpu
}//namespace exblas
} //namespace dg
#endif // EXSUM_FPE_HPP_
