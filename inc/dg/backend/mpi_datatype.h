#pragma once

#include <mpi.h>
#include <complex>
#include <type_traits>
#include <thrust/complex.h>
#include "predicate.h"

namespace dg
{
///@cond
// The MPI type corresponding to a C++ type (unsupported types do not compile)
template<class value_type>
inline MPI_Datatype getMPIDataType()
{
    using T = value_type;
    if constexpr( std::is_same_v<T, char>) return MPI_CHAR;
    else if constexpr( std::is_same_v<T, signed char>) return MPI_SIGNED_CHAR;
    else if constexpr( std::is_same_v<T, unsigned char>) return MPI_UNSIGNED_CHAR;
    else if constexpr( std::is_same_v<T, wchar_t>) return MPI_WCHAR;
    else if constexpr( std::is_same_v<T, short>) return MPI_SHORT;
    else if constexpr( std::is_same_v<T, unsigned short>) return MPI_UNSIGNED_SHORT;
    else if constexpr( std::is_same_v<T, int>) return MPI_INT;
    else if constexpr( std::is_same_v<T, unsigned>) return MPI_UNSIGNED;
    else if constexpr( std::is_same_v<T, long>) return MPI_LONG;
    else if constexpr( std::is_same_v<T, unsigned long>) return MPI_UNSIGNED_LONG;
    else if constexpr( std::is_same_v<T, long long>) return MPI_LONG_LONG;
    else if constexpr( std::is_same_v<T, unsigned long long>) return MPI_UNSIGNED_LONG_LONG;
    else if constexpr( std::is_same_v<T, float>) return MPI_FLOAT;
    else if constexpr( std::is_same_v<T, double>) return MPI_DOUBLE;
    else if constexpr( std::is_same_v<T, long double>) return MPI_LONG_DOUBLE;
    else if constexpr( std::is_same_v<T, bool>) return MPI_C_BOOL;
    else if constexpr( std::is_same_v<T, std::complex<float>> or
                       std::is_same_v<T, thrust::complex<float>>) return MPI_C_COMPLEX;
    else if constexpr( std::is_same_v<T, std::complex<double>> or
                       std::is_same_v<T, thrust::complex<double>>) return MPI_C_DOUBLE_COMPLEX;
    else if constexpr( std::is_same_v<T, std::complex<long double>> or
                       std::is_same_v<T, thrust::complex<long double>>) return MPI_C_LONG_DOUBLE_COMPLEX;
    else
        static_assert( detail::always_false_v<T>, "Type not supported by MPI!");
}
///@endcond
}//namespace dg
