/*
 * %%%%%%%%%%%%%%%%%%%%%%%Original development%%%%%%%%%%%%%%%%%%%%%%%%%
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 *  All rights reserved.
 * %%%%%%%%%%%%%%%%%%%%%%%Modifications and further additions%%%%%%%%%%
 *  Matthias Wiesenberger, 2017, within FELTOR and EXBLAS licenses
 */
/**
 *  @file accumulate.h
 *  @brief Primitives for accumulation into superaccumulator
 *
 *  @authors
 *    Developers : \n
 *        Roman Iakymchuk  -- roman.iakymchuk@lip6.fr \n
 *        Sylvain Collange -- sylvain.collange@inria.fr \n
 *        Matthias Wiesenberger -- mattwi@fysik.dtu.dk
 */
#pragma once
#include <array>
#include "config.h"
#include "mylibm.hpp"
//this file has a direct correspondance to gpu code accumulate.cuh

namespace dg
{
namespace exblas {
namespace cpu {
///@cond
template<class T>
inline T get_element( T x, int){
	return x;
}
template<class T>
inline T get_element( T* x, int i){
	return *(x+i);
}
// The k-th real component of a value that the exact sums accumulate in the
// k-th superaccumulator (see ValueTraits)
inline double get_result( double x, unsigned){
    return x;
}
inline double get_result( float x, unsigned){
    return double(x);
}
template<class real>
inline double get_result( std::complex<real> x, unsigned k){
    return k == 0 ? double(x.real()) : double(x.imag());
}
template<class real>
inline double get_result( thrust::complex<real> x, unsigned k){
    return k == 0 ? double(x.real()) : double(x.imag());
}
////////////////////////////////////////////////////////////////////////////////
// Auxiliary functions
////////////////////////////////////////////////////////////////////////////////

// Knuth 2Sum.
template<typename T>
inline std::enable_if_t<!std::is_integral_v<T>,T> KnuthTwoSum(T a, T b, T & s)
{
    T r = a + b;
    T z = r - a;
    s = (a - (r - z)) + (b - z);
    return r;
}

////////////////////////////////////////////////////////////////////////////////
// Main computation pass: compute partial superaccs
////////////////////////////////////////////////////////////////////////////////
inline void AccumulateWord( int64_t *accumulator, int i, int64_t x) {
    // With atomic accumulator updates
    // accumulation and carry propagation can happen in any order,
    // as long as addition is atomic
    // only constraint is: never forget an overflow bit
    unsigned char overflow;
    int64_t carry = x;
    int64_t carrybit;
    int64_t oldword = cpu::xadd(accumulator[i], x, overflow);
    // To propagate over- or underflow
    while(EXBLAS_UNLIKELY(overflow)) {
        // Carry or borrow
        // oldword has sign S
        // x has sign S
        // accumulator[i] has sign !S (just after update)
        // carry has sign !S
        // carrybit has sign S
        // unsigned addition and shift: signed overflow and left shift of
        // negative numbers are undefined
        carry = (int64_t)((uint64_t)oldword + (uint64_t)carry) >> DIGITS;    // Arithmetic shift
        bool s = oldword > 0;
        carrybit = (s ? 1ll << KRX : (unsigned long long)(-1ll) << KRX); //MW left shift of negative number is undefined so convert to unsigned

        // Cancel carry-save bits
        cpu::xadd(accumulator[i], (int64_t)( 0 - ((uint64_t)carry << DIGITS)), overflow);
        if(EXBLAS_TSAFE && EXBLAS_UNLIKELY(s ^ overflow)) {
            // (Another) overflow of sign S
            carrybit *= 2;
        }
        carry += carrybit;

        ++i;
        if (i >= BIN_COUNT){
            return;
        }
        oldword = cpu::xadd(accumulator[i], carry, overflow);
    }
}
///@endcond

/**
* @brief Accumulate a double to the superaccumulator
*
* @param accumulator a pointer to at least \c BIN_COUNT 64 bit integers on the CPU (representing the superaccumulator)
* @param x the double to add to the superaccumulator
*/
inline void Accumulate( int64_t* accumulator, double x) {
    if (x == 0)
        return;


    int e = cpu::exponent(x);
    int exp_word = e / DIGITS;  // Word containing MSbit (upper bound)
    int iup = exp_word + F_WORDS;

    double xscaled = cpu::myldexp(x, -DIGITS * exp_word);

    int i;
    for (i = iup; i>=0 && xscaled != 0; --i) { //MW: i>=0 protects against NaN
        double xrounded = cpu::myrint(xscaled);
        int64_t xint = cpu::myllrint(xscaled);
        AccumulateWord(accumulator, i, xint);

        xscaled -= xrounded;
        xscaled *= DELTASCALE;
    }
}
#ifndef EXBLAS_WITHOUT_VCL
inline void Accumulate( int64_t* accumulator, vcl::Vec8d x) {
    double v[8];
    x.store(v);

#if INSTRSET >= 7
    _mm256_zeroupper();
#endif
    for(unsigned int j = 0; j != 8; ++j) {
        exblas::cpu::Accumulate(accumulator, v[j]);
    }
}
#endif //EXBLAS_WITHOUT_VCL
////////////////////////////////////////////////////////////////////////////////
// Normalize functions
////////////////////////////////////////////////////////////////////////////////
// Returns sign
// Does not really normalize! MW: what does that mean?
//
/**
* @brief Normalize a superaccumulator
*
* @param accumulator a pointer to at least \c BIN_COUNT 64 bit integers on the CPU (representing the superaccumulator)
* @param imin the first index in the accumulator
* @param imax the last index in the accumulator
*
* @return  carry in bit (sign)
*/
inline bool Normalize( int64_t *accumulator, int& imin, int& imax) {
    // unsigned left shift: left shift of negative numbers is undefined
    int64_t carry_in = accumulator[imin] >> DIGITS;
    accumulator[imin] -= (int64_t)((uint64_t)carry_in << DIGITS);
    int i;
    // Sign-extend all the way
    for (i = imin + 1; i < BIN_COUNT; ++i) {
        accumulator[i] += carry_in;
        int64_t carry_out = accumulator[i] >> DIGITS;    // Arithmetic shift
        accumulator[i] -= (int64_t)((uint64_t)carry_out << DIGITS);
        carry_in = carry_out;
    }
    imax = i - 1;

    // Do not cancel the last carry to avoid losing information
    accumulator[imax] += (int64_t)((uint64_t)carry_in << DIGITS);

    return carry_in < 0;
}

/**
* @brief Add superaccumulators: \c acc+=x
*
* Both are normalized first, so that the sum of their words cannot overflow
* @param acc pointer to \c num_superacc*BIN_COUNT 64 bit integers on the CPU
* (\c num_superacc superaccumulators stored one after the other)
* @param x pointer to \c num_superacc*BIN_COUNT 64 bit integers on the CPU
* (overwritten with its normalized form)
* @param num_superacc number of superaccumulators in \c acc and \c x
*/
inline void AddSuperaccs( int64_t* acc, int64_t* x, unsigned num_superacc)
{
    for( unsigned k=0; k<num_superacc; k++)
    {
        int imin = IMIN, imax = IMAX;
        Normalize( acc + k*BIN_COUNT, imin, imax);
        imin = IMIN, imax = IMAX;
        Normalize( x + k*BIN_COUNT, imin, imax);
        for( int i=IMIN; i<=IMAX; i++)
            acc[k*BIN_COUNT+i] += x[k*BIN_COUNT+i];
    }
}

////////////////////////////////////////////////////////////////////////////////
// Rounding functions
////////////////////////////////////////////////////////////////////////////////


/**
* @brief Convert a superaccumulator to the nearest double precision number (CPU version)
*
* @tparam T \c double (the default), \c float or a complex number of these
* (\c std::complex or \c thrust::complex). A \c float is the \c double
* result converted to \c float. A complex number is made of the results of
* two superaccumulators, the one of the real part first, as \c exsum_cpu, \c
* exsum_omp and \c exsum_gpu write them.
* @param accumulator a pointer to at least \c BIN_COUNT 64 bit integers on
* the CPU (representing the superaccumulator), or \c 2*BIN_COUNT for complex
* \c T (normalized on output)
* @return the double precision number nearest to the superaccumulator,
* converted to \c T
*/
template<class T = double>
inline T Round( int64_t * accumulator) {
    if constexpr( ValueTraits<T>::K == 2)
    {
        using real = typename T::value_type;
        return T( real( Round<double>( accumulator)),
                  real( Round<double>( accumulator + BIN_COUNT)));
    }
    else if constexpr( !std::is_same_v<T, double>)
        return T( Round<double>( accumulator));
    else
    {
        int imin = IMIN;
        int imax = IMAX;
        bool negative = Normalize(accumulator, imin, imax);

        // Find leading word
        int i;
        // Skip zeroes
        for(i = imax; i >= imin && accumulator[i] == 0; --i) {
            // MW: note that i >= imin has to come *before* accumulator[i]
            // else it is possible that accumulator[-1] is accessed
        }
        if (negative) {
            // Skip ones
            for(; i >= imin && (accumulator[i] & ((1ll << DIGITS) - 1)) == ((1ll << DIGITS) - 1); --i) {
            }
        }
        if (i < 0) {
            return 0.0;
        }

        int64_t hiword = negative ? ((1ll << DIGITS) - 1) - accumulator[i] : accumulator[i];
        double rounded = (double)hiword;
        double hi = ldexp(rounded, (i - F_WORDS) * DIGITS);
        if (i == 0) {
            return negative ? -hi : hi;  // Correct rounding achieved
        }
        hiword -= std::llrint(rounded);
        double mid = ldexp((double) hiword, (i - F_WORDS) * DIGITS);

        // Compute sticky
        int64_t sticky = 0;
        for (int j = imin; j != i - 1; ++j) {
            sticky |= negative ? ((1ll << DIGITS) - accumulator[j]) : accumulator[j];
        }

        int64_t loword = negative ? ((1ll << DIGITS) - accumulator[i - 1]) : accumulator[i - 1];
        loword |= !!sticky;
        double lo = ldexp((double) loword, (i - 1 - F_WORDS) * DIGITS);

        // Now add3(hi, mid, lo)
        // No overlap, we have already normalized
        if (mid != 0) {
            lo = cpu::OddRoundSumNonnegative(mid, lo);
        }
        // Final rounding
        hi = hi + lo;
        return negative ? -hi : hi;
    }
}

}//namespace cpu
} //namespace exblas
} //namespace dg
