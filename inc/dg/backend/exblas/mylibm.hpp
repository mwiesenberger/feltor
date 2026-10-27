/*
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 */

/**
 *  \file mylibm.hpp
 *  \brief Provides a set of auxiliary functions to superaccumulation.
 *         For internal use
 *
 *  \authors
 *    Developers : \n
 *        Roman Iakymchuk  -- roman.iakymchuk@lip6.fr \n
 *        Sylvain Collange -- sylvain.collange@inria.fr \n
 *        Matthias Wiesenberger -- mattwi@fysik.dtu.dk
 */

#ifndef MYLIBM_HPP_INCLUDED
#define MYLIBM_HPP_INCLUDED

#include <cstring> // std::memcpy
#include "config.h"

namespace dg
{
namespace exblas{
namespace cpu{

inline int64_t myllrint(double x) {
#if not defined EXBLAS_WITHOUT_VCL && not defined( _MSC_VER)
    return _mm_cvtsd_si64(_mm_set_sd(x));
#else
    return std::llrint(x);
#endif
}

inline double myrint(double x)
{
#ifndef EXBLAS_WITHOUT_VCL
#if defined( __GNUG__) || defined( _MSC_VER)
    // Workaround gcc bug 51033
    //_mm_round_sd(_mm_undefined_pd(), _mm_set_sd(x));
    //__m128d undefined;
    //r.v = _mm_round_sd(_mm_setzero_pd(), _mm_set_sd(x), _MM_FROUND_TO_NEAREST_INT);
    //r.v = _mm_round_sd(undefined, _mm_set_sd(x), _MM_FROUND_TO_NEAREST_INT);
    return _mm_cvtsd_f64( _mm_round_pd(_mm_set_sd(x), _MM_FROUND_TO_NEAREST_INT));
#else
    double r;
    //asm("roundsd $0, %1, %0" : "=x" (r) : "x" (x));
    asm(EXBLAS_ASM_BEGIN "roundsd %0, %1, 0" EXBLAS_ASM_END : "=x" (r) : "x" (x));
    return r;
#endif
#else
    return std::rint(x);
#endif//EXBLAS_WITHOUT_VCL
}

// Reinterpret the bits of a double; unlike type punning with a union,
// std::memcpy is well-defined in C++ and compiles to a register move
inline uint64_t double_as_uint64(double x)
{
    uint64_t i;
    std::memcpy( &i, &x, sizeof(i));
    return i;
}
inline double uint64_as_double(uint64_t i)
{
    double x;
    std::memcpy( &x, &i, sizeof(x));
    return x;
}

inline int exponent(double x)
{
    // simpler frexp
    uint64_t e = ((double_as_uint64(x) >> 52) & 0x7ff) - 0x3ff;
    return e;
}

//MW: if x is a NaN this thing still returns a number??
inline double myldexp(double x, int e)
{
    // Scale x by e
    uint64_t i = double_as_uint64(x);
    // Subnormal numbers have no implicit leading bit, so adding to the
    // exponent bits does not scale them
    if( EXBLAS_UNLIKELY( (i & 0x7ff0000000000000ull) == 0))
        return std::ldexp( x, e);
    i += (uint64_t)e << 52;
    return uint64_as_double(i);
}

// Assumptions: th>tl>=0, no overlap between th and tl
inline double OddRoundSumNonnegative(double th, double tl)
{
    // Adapted from:
    // Sylvie Boldo, and Guillaume Melquiond. "Emulation of a FMA and correctly rounded sums: proved algorithms using rounding to odd." IEEE Transactions on Computers, 57, no. 4 (2008): 462-471.
    uint64_t thdb = double_as_uint64( th + tl);
    // - if the mantissa of th is odd, there is nothing to do
    // - otherwise, round up as both tl and th are positive
    // in both cases, this means setting the msb to 1 when tl>0
    thdb |= (tl != 0.0);
    return uint64_as_double( thdb);
}

#ifdef THREADSAFE //MW: we don't define this anywhere?
#define EXBLAS_TSAFE 1
#define EXBLAS_LOCK_PREFIX "lock "
#else
#define EXBLAS_TSAFE 0 //MW: so I guess it's always 0 i.e. false
#define EXBLAS_LOCK_PREFIX
#endif

// signedcarry in {-1, 0, 1}
inline int64_t xadd(int64_t & memref, int64_t x, unsigned char & of)
{

//msvc doesn't allow inline assembler code
//If we don't have VCL, then sometimes the assembler code also makes problems
#if (defined (EXBLAS_WITHOUT_VCL) || defined(_MSC_VER)) && !EXBLAS_TSAFE
//manually compute non-atomic load-ADDC-store
	int64_t y = memref;
	// unsigned addition: signed overflow is undefined
	memref = (int64_t)((uint64_t)y + (uint64_t)x);
	int64_t x63 = (x >> 63) & 1;
	int64_t y63 = (y >> 63) & 1;
	int64_t r63 = (memref >> 63) & 1;
	int64_t c62 = r63 ^ x63 ^ y63;
	int64_t c63 = (x63 & y63) | (c62 & (x63 | y63));
	of = c63 ^ c62;
	return y;
#else
    // OF and SF  -> carry=1
    // OF and !SF -> carry=-1
    // !OF        -> carry=0
    int64_t oldword = x;
#ifdef ATT_SYNTAX
    asm volatile (EXBLAS_LOCK_PREFIX"xaddq %1, %0\n"
        "setob %2"
     : "+m" (memref), "+r" (oldword), "=q" (of) : : "cc", "memory");
#else
    asm volatile (EXBLAS_LOCK_PREFIX"xadd %1, %0\n"
        "seto %2"
     : "+m" (memref), "+r" (oldword), "=q" (of) : : "cc", "memory");
#endif //ATT_SYNTAX
    return oldword;
#endif //_MSC_VER
}

#ifndef EXBLAS_WITHOUT_VCL
inline bool horizontal_or(vcl::Vec8d const & a) {
    //return _mm512_movemask_pd(a) != 0;
    vcl::Vec8db p = a != 0;
    return vcl::horizontal_or( p);
    //return !_mm512_testz_pd(p, p);
}
#else
inline bool horizontal_or( const double & a){
    return a!= 0;
}
#endif//EXBLAS_WITHOUT_VCL


}//namespace cpu
}//namespace exblas
} // namespace dg

#endif
