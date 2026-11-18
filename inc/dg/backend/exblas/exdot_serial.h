/*
 * %%%%%%%%%%%%%%%%%%%%%%%Original development%%%%%%%%%%%%%%%%%%%%%%%%%
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 * %%%%%%%%%%%%%%%%%%%%%%%Modifications and further additions%%%%%%%%%%
 *  Matthias Wiesenberger, 2017, within FELTOR and EXBLAS licenses
 */

/**
 *  @file exdot_serial.h
 *  @brief Serial version of exdot
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

#include "accumulate.h"
#include "ExSUM.FPE.hpp"

namespace dg
{
namespace exblas{

/*! @brief Utility union to display all bits of a double (using <a href="https://en.wikipedia.org/wiki/Type_punning">type-punning</a>)
@warning Reading the member that was not written last is undefined behaviour
in C++ (unlike in C). GCC, Clang and MSVC support it as an extension, but for
well-defined code copy the bits with \c std::memcpy (or \c std::bit_cast in
C++20) instead

@code
double result; // = ...
udouble res;
res.d = result;
std::cout << "Result as double "<<res.d<<"  as integer "<<res.i<<std::endl;
@endcode
*/
union udouble{
    double d; //!< a double
    int64_t i; //!< a 64 bit integer
};

/*! @brief Utility union to display all bits of a float (using <a href="https://en.wikipedia.org/wiki/Type_punning">type-punning</a>)
@warning Reading the member that was not written last is undefined behaviour
in C++ (unlike in C). GCC, Clang and MSVC support it as an extension, but for
well-defined code copy the bits with \c std::memcpy (or \c std::bit_cast in
C++20) instead

@code
float result; // = ...
ufloat res;
res.f = result;
std::cout << "Result as float "<<res.f<<"  as integer "<<res.i<<std::endl;
@endcode
*/
union ufloat{
    float f; //!< a float
    int32_t i; //!< a 32 bit integer
};

///@cond
namespace cpu{

// Exact sum of f( x_0[i], x_1[i], ...) into the superaccumulator acc
// With VCL the FPE works on 8 values at once. This hand vectorisation
// matters: the compiler does not vectorise the scalar loop, since every value
// goes through a chain of TwoSums into the same FPE (a loop-carried
// dependency) and a data dependent flush. Vectorising needs 8 independent FPEs,
// which changes the order of additions; only the exactness of the sum makes
// that allowed. Measured on an i5-14400F, serial x*y: 30-36 GB/s with VCL vs
// 4-5 GB/s for the scalar loop.
// f is called on the elements and its 8 values are collected into a
// vcl::Vec8d (instead of calling f on vcl::Vec8d), so it computes in the types
// of the elements for every functor, like the GPU version does. This is not
// slower than computing a product on vcl::Vec8d
template<size_t NBFPE, class Functor, class ...PointerOrValues>
void ExSUMFPE_cpu(int N, int64_t* acc, bool* error, Functor f, PointerOrValues ...xs) {
    using value_type = return_type<Functor, PointerOrValues...>;
    constexpr unsigned K = ValueTraits<value_type>::K; // one superaccumulator per component
#ifndef EXBLAS_WITHOUT_VCL
    using Cache = FPExpansionVect<vcl::Vec8d, NBFPE, FPExpansionTraits<true> >;
    auto cache = make_caches<Cache>( [acc]( size_t c){ return acc + c*BIN_COUNT;},
        std::make_index_sequence<K>());
    int r = (( int64_t(N) ) & ~7ul);
    for(int i = 0; i < r; i+=8) {
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
            if( !vcl::horizontal_and( finite) ) *error = true;
            cache[c].Accumulate(x);
        }
    }
    if( r != N) {
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
            if( !vcl::horizontal_and( finite) ) *error = true;
            cache[c].Accumulate(x);
        }
    }
#else// EXBLAS_WITHOUT_VCL
    using Cache = FPExpansionVect<double, NBFPE, FPExpansionTraits<true> >;
    auto cache = make_caches<Cache>( [acc]( size_t c){ return acc + c*BIN_COUNT;},
        std::make_index_sequence<K>());
    for(int i = 0; i < N; i++) {
        value_type v = f( get_element(xs,i)...);
        for( unsigned c=0; c<K; c++)
        {
            double x = ValueTraits<value_type>::get( v, c);
            if( !std::isfinite(x) ) *error = true;
            cache[c].Accumulate(x);
        }
    }
#endif// EXBLAS_WITHOUT_VCL
    for( unsigned c=0; c<K; c++)
        cache[c].Flush();
}
}//namespace cpu
///@endcond

/*!@class hide_hostacc
 * @param h_superacc pointer to an array of 64 bit integers (the
 * superaccumulator) in **host memory** with size at least \c exblas::BIN_COUNT
 * (39) for double and float values, or \c 2*exblas::BIN_COUNT for complex
 * values of \c f, where the superaccumulator of the real part comes first
 * (contents are overwritten, the function does not allocate memory i.e. the
 * memory needs to be allocated **before** calling the function)
 * @param status 0 indicates success, 1 indicates that a value of \c f was
 * NaN or Inf (which may be caused by NaN or Inf inputs to \c f)
 */
/*!@class hide_gpu_static_buffer
 * @attention The function uses a static buffer in device memory (and possibly
 * in pinned host memory), which is allocated on the device that is current at
 * the first call. Therefore
 * - it must not be called concurrently from several host threads (e.g.
 *   OpenMP threads), even if they use the same device, because all calls
 *   share the buffer, and
 * - it must only be called while that same device is current (so e.g. one
 *   host thread per GPU is not possible)
 *
 * One MPI process per GPU is possible, since every process has its own buffer
 */

/*!@class hide_exsum
 *
 * Accumulate the exact sum \f[ \sum_{i=0}^{N-1} f( x_{0i}, x_{1i}, \dots) \f] into a superaccumulator.
 * For double and float values of \c f the superaccumulator is an array of \c
 * exblas::BIN_COUNT (39) 64 bit integers that represents a large fixed point
 * number such that the summation is computed with virtually infinite precision
 * and is thus bitwise reproducible even in a parallel environment. For complex
 * \c f two separate superaccumulators are used, one for the real and one for
 * the imaginary parts of \c f.
 * @note Any superaccumulator can be converted to a double precision number
 * using the Round function: \ref dg::exblas::cpu::Round
 * @attention the values of \c f are **not** computed with infinite precision
 * only the sum is (this does not break the reproducibility)
 * @note The algorithm is described in the paper
 <a href = "https://hal.archives-ouvertes.fr/hal-00949355v3">Sylvain Collange, David Defour, Stef Graillat, Roman Iakymchuk. "Numerical Reproducibility for the Parallel Reduction on Multi- and Many-Core Architectures", 2015. </a>
 * @tparam NBFPE size of the floating point expansion that caches the values
 * before they are added to the superaccumulator (typically 3 on the GPU and 8 on CPU code)
 * @tparam Functor signature <tt> return_type operator()( value_type_x0, value_type_x1, ...) </tt>
 * where \c return_type is \c float, \c double or a complex number of these (\c
 * std::complex or \c thrust::complex). float values in return_type are
 * converted to \c double when added to the superaccumulator.  \c f must be
 * callable on the device in use.
 * @tparam PointerOrValues each must be a value \c T or a pointer <tt> T* </tt>
 * or <tt> const T* </tt> with T convertible to the respective argument type
 * of \c f. If it is a pointer type, then we iterate through the pointed data
 * from 0 to \c size, else we consider the value constant in every iteration.
 * Pointers must point to device memory for \c exsum_gpu and to host memory
 * for \c exsum_cpu and \c exsum_omp.
 * @param size size N of the arrays to sum
 * @param f the functor to evaluate
 * @param xs the arrays
 */

///@brief Serial version of exact sum
///@copydoc hide_exsum
///@copydoc hide_hostacc
template<size_t NBFPE=8, class Functor, class ...PointerOrValues>
void exsum_cpu(unsigned size, int64_t* h_superacc, int* status, Functor f, PointerOrValues ...xs){
    constexpr unsigned K = ValueTraits<return_type<Functor, PointerOrValues...>>::K;
    for( unsigned i=0; i<K*exblas::BIN_COUNT; i++)
        h_superacc[i] = 0;
    bool error = false;
    cpu::ExSUMFPE_cpu<NBFPE>( (int)size, h_superacc, &error, f, xs...);
    *status = 0;
    if( error ) *status = 1;
}

}//namespace exblas
} //namespace dg
