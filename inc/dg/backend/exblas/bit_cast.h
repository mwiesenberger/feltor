#pragma once
/**
 *  @file bit_cast.h
 *  @brief The C++20 function bit_cast for C++17
 */
#include <cstring> // std::memcpy
#include <type_traits>

namespace dg
{
namespace exblas
{
/*! @brief Reinterpret the bits of \c src as a \c To, e.g. to compare results bitwise
 *
 * The same as \c std::bit_cast of C++20 (except that it is not \c
 * constexpr), so that with C++20 \c dg::exblas::bit_cast can be replaced by \c
 * std::bit_cast. Unlike reading \c udouble::i after writing \c udouble::d it
 * copies the bits with \c std::memcpy, which is well-defined C++
@code
double result; // = ...
std::cout << "Result as double "<<result<<"  as integer "<<dg::exblas::bit_cast<int64_t>( result)<<std::endl;
@endcode
 * @tparam To the destination type, trivially copyable and trivially constructible
 * @tparam From the source type, trivially copyable and of the size of \c To
 * @param src the object whose bits are copied
 * @return an object of type \c To with the bits of \c src
 */
template<class To, class From>
std::enable_if_t<
    sizeof(To) == sizeof(From) &&
    std::is_trivially_copyable_v<From> &&
    std::is_trivially_copyable_v<To>,
    To>
bit_cast(const From& src) noexcept
{
    static_assert(std::is_trivially_constructible_v<To>,
        "This implementation additionally requires "
        "destination type to be trivially constructible");

    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

}//namespace exblas
} //namespace dg
