#pragma once

#include <variant>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <sstream>
#include <iomanip> // std::put_time
#include <ctime>  // std::localtime
#include <netcdf.h>
#include "nc_error.h"

namespace dg
{
namespace file
{

/*! @brief Utility type to simplify dealing with heterogeneous attribute types
 *
 *  @note Unfortunately, user defined types exist so not every attribute can be
 *  an nc_att_t
 *  @ingroup netcdf
*/
using nc_att_t = std::variant<int, unsigned, float, double, bool, std::string,
      std::vector<int>, std::vector<unsigned>, std::vector<float>,
      std::vector<double>, std::vector<bool>>;

/*! @brief Generate one line entry for the history global attribute
 *
 * This will generate a string containing a whitespace seperated list of
 *  -# the current *UTC* time in ISO 8601 format: <tt>%Y-%m-%dT%H:%M:%SZ</tt>
 *  -# all given argv (whitespace separated)
 *  -# A newline
 *  .
 * @snippet easy_atts_t.cpp timestamp
 * @param argc from main function
 * @param argv from main function
 * @return string containing current time followed by the parameters with which the program was invoked
 * @sa See history in <a href="https://docs.unidata.ucar.edu/netcdf-c/current/attribute_conventions.html">Attribute Convenctions</a>
 * @sa <a href="https://en.cppreference.com/w/cpp/io/manip/put_time">std::put_time</a>
 * @ingroup netcdf
 */
inline std::string timestamp( int argc, char* argv[])
{
    // TODO: once we're using C++20 this should be implemented using chrono
    // Get UTC time
    auto ttt = std::time(nullptr);
    std::ostringstream oss;
    // time string  + program-name + args
    oss << std::put_time(std::gmtime(&ttt), "%FT%TZ");
    for( int i=0; i<argc; i++) oss << " "<<argv[i];
    oss << std::endl;
    return oss.str();
}

/*! @brief Version compile time constants available as a map
 *
 * Is intended to be used as NetCDF file attributes
 * @snippet easy_atts_t.cpp timestamp
 * @note We use underscore instead of "git-hash"  so that python netcdf can
 *   more easily read the attribute
 *
 * The entries in the map are filled only if the corresponding MACROs are defined at compile time.
 * Use <tt>-DVERSION_FLAGS</tt> during compilation to define all otherwise it remains empty.
 * This is the corresponding entry in \c feltor/config/version.mk
 * @code{.sh}
 * GIT_HASH=$(git rev-parse HEAD)
 * COMPILE_TIME=$(date -u +'%Y-%m-%dT%H:%M:%SZ')
 * GIT_BRANCH=$(git branch --show-current)
 * @endcode
 * @sa This approach follows
 <a href="https://stackoverflow.com/questions/44038428/include-git-commit-hash-and-or-branch-name-in-c-c-source/44038455#44038455">stackoverflow</a>.
 *  @ingroup netcdf
 */
inline const std::map<std::string, std::string> version_flags =
{
#ifdef GIT_HASH
    {"git_hash", GIT_HASH},
#endif // GIT_HASH
#ifdef GIT_BRANCH
    {"git_branch", GIT_BRANCH},
#endif // GIT_BRANCH
#ifdef COMPILE_TIME
    {"compile_time", COMPILE_TIME},
#endif // COMPILE_TIME
};

///@cond
namespace detail
{
template<class T>
inline constexpr bool always_false_v = false;

// The NetCDF type corresponding to a C++ type (unsupported types do not compile)
template<class value_type>
constexpr nc_type getNCDataType()
{
    using T = value_type;
    if constexpr( std::is_same_v<T, double>) return NC_DOUBLE;
    else if constexpr( std::is_same_v<T, float>) return NC_FLOAT;
    else if constexpr( std::is_same_v<T, char>) return NC_CHAR;
    else if constexpr( std::is_same_v<T, signed char>) return NC_BYTE;
    else if constexpr( std::is_same_v<T, unsigned char>) return NC_UBYTE;
    else if constexpr( std::is_same_v<T, short>) return NC_SHORT;
    else if constexpr( std::is_same_v<T, unsigned short>) return NC_USHORT;
    else if constexpr( std::is_same_v<T, int>) return NC_INT;
    else if constexpr( std::is_same_v<T, unsigned>) return NC_UINT;
    // long is 64 bit on Linux and macOS but 32 bit on Windows
    else if constexpr( std::is_same_v<T, long>) return sizeof(long) == 8 ? NC_INT64 : NC_INT;
    else if constexpr( std::is_same_v<T, long long>) return NC_INT64;
    else if constexpr( std::is_same_v<T, unsigned long long>) return NC_UINT64;
    else if constexpr( std::is_same_v<T, bool>) return NC_BYTE;
    else if constexpr( std::is_same_v<T, std::string> or
                       std::is_same_v<T, const char*>) return NC_STRING;
    else
        static_assert( always_false_v<T>, "Type not supported by NetCDF!");
}

// C types for which NetCDF has typed attribute functions nc_put_att_<type> and
// nc_get_att_<type> (strings are handled separately)
template<class T>
inline constexpr bool is_nc_att_type_v =
    std::is_same_v<T, char> or std::is_same_v<T, signed char> or
    std::is_same_v<T, unsigned char> or std::is_same_v<T, short> or
    std::is_same_v<T, unsigned short> or std::is_same_v<T, int> or
    std::is_same_v<T, unsigned> or std::is_same_v<T, long> or
    std::is_same_v<T, long long> or std::is_same_v<T, unsigned long long> or
    std::is_same_v<T, float> or std::is_same_v<T, double>;

// Call the NetCDF function matching the C type T. NetCDF converts the values
// to the attribute type xtype (and returns e.g. NC_ERANGE if a value is not
// representable). We do not add any checks of our own.
template<class T>
int put_att_T( int ncid, int varid, const char* name, nc_type xtype,
    size_t size, const T* data)
{
    static_assert( is_nc_att_type_v<T>);
    if constexpr( std::is_same_v<T, char>) // always NC_CHAR
        return nc_put_att_text( ncid, varid, name, size, data);
    else if constexpr( std::is_same_v<T, signed char>)
        return nc_put_att_schar( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, unsigned char>)
        return nc_put_att_ubyte( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, short>)
        return nc_put_att_short( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, unsigned short>)
        return nc_put_att_ushort( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, int>)
        return nc_put_att_int( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, unsigned>)
        return nc_put_att_uint( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, long>)
        return nc_put_att_long( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, long long>)
        return nc_put_att_longlong( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, unsigned long long>)
        return nc_put_att_ulonglong( ncid, varid, name, xtype, size, data);
    else if constexpr( std::is_same_v<T, float>)
        return nc_put_att_float( ncid, varid, name, xtype, size, data);
    else
        return nc_put_att_double( ncid, varid, name, xtype, size, data);
}

// Call the NetCDF function matching the C type T. NetCDF converts the values
// from the attribute type (and returns e.g. NC_ERANGE if a value is not
// representable or NC_ECHAR for text <-> number). No checks of our own.
template<class T>
int get_att_T( int ncid, int varid, const char* name, T* data)
{
    static_assert( is_nc_att_type_v<T>);
    if constexpr( std::is_same_v<T, char>)
        return nc_get_att_text( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, signed char>)
        return nc_get_att_schar( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, unsigned char>)
        return nc_get_att_ubyte( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, short>)
        return nc_get_att_short( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, unsigned short>)
        return nc_get_att_ushort( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, int>)
        return nc_get_att_int( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, unsigned>)
        return nc_get_att_uint( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, long>)
        return nc_get_att_long( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, long long>)
        return nc_get_att_longlong( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, unsigned long long>)
        return nc_get_att_ulonglong( ncid, varid, name, data);
    else if constexpr( std::is_same_v<T, float>)
        return nc_get_att_float( ncid, varid, name, data);
    else
        return nc_get_att_double( ncid, varid, name, data);
}

// S allows both std::string and const char* to be used
template<class S, class T>
void put_att( int ncid, int varid, const std::tuple<S, nc_type,
        std::vector<T>>& att)
{
    // This will convert const char* to std::string
    std::string name_string = std::get<0>(att);
    auto name = name_string.c_str();
    nc_type xtype = std::get<1>(att);
    const std::vector<T>& data = std::get<2>(att);
    unsigned size = data.size();
    NC_Error_Handle err;
    if constexpr( is_nc_att_type_v<T>)
    {
        err = put_att_T( ncid, varid, name, xtype, size, data.data());
    }
    else if constexpr( std::is_same_v<T, std::string> or
                       std::is_same_v<T, const char*>)
    {
        if( xtype == NC_STRING) // any number of strings
        {
            std::vector<const char*> strs( size);
            for( unsigned u=0; u<size; u++)
            {
                if constexpr( std::is_same_v<T, std::string>)
                    strs[u] = data[u].c_str();
                else
                    strs[u] = data[u];
            }
            err = nc_put_att_string( ncid, varid, name, size, strs.data());
        }
        else if( xtype == NC_CHAR) // exactly one string
        {
            if( size != 1)
                throw std::runtime_error( "Cannot write a string array attribute "
                    +name_string+" as NC_CHAR (use NC_STRING)");
            std::string tmp = data[0];
            err = nc_put_att_text( ncid, varid, name, tmp.size(), tmp.c_str());
        }
        else
            throw std::runtime_error( "String attribute "+name_string
                +" can only be written as NC_CHAR or NC_STRING");
    }
    else if constexpr( std::is_same_v<T, bool>)
    {
        // std::vector<bool> is not necessarily contiguous
        std::vector<signed char> dataB(size);
        for( unsigned i=0; i<size; i++)
            dataB[i] = data[i];
        err = nc_put_att_schar( ncid, varid, name, xtype, size, dataB.data());
    }
    else // user defined type (e.g. compound), written without conversion
    {
        size_t type_size;
        err = nc_inq_type( ncid, xtype, NULL, &type_size);
        if( type_size != sizeof(T))
            throw std::runtime_error( "Size of given type does not match "
                "NetCDF type of attribute "+name_string);
        err = nc_put_att( ncid, varid, name, xtype, size, data.data());
    }
}

template<class S, class T> // T cannot be nc_att_t
void put_att( int ncid, int varid, std::tuple<S, nc_type, T> att)
{
    put_att( ncid, varid, std::make_tuple( std::get<0>(att), std::get<1>(att),
                std::vector<T>( 1, std::get<2>(att)) ));
}

// Type of an attribute if none is given explicitly: a single string is
// written as NC_CHAR (the classic text attribute), several as NC_STRING
template<class T>
nc_type default_att_type( size_t size)
{
    if constexpr( std::is_same_v<T, std::string> or
                  std::is_same_v<T, const char*>)
        return size == 1 ? NC_CHAR : NC_STRING;
    else
        return getNCDataType<T>();
}

// Variants for normal types
template<class S, class T>
void put_att( int ncid, int varid, const std::pair<S, T>& att)
{
    put_att( ncid, varid, std::make_tuple( att.first,
                detail::default_att_type<T>(1), std::vector<T>(1,att.second)));
}

template<class S, class T>
void put_att( int ncid, int varid, const std::pair<S, std::vector<T>>& att)
{
    put_att( ncid, varid, std::make_tuple( att.first,
                detail::default_att_type<T>( att.second.size()), att.second));
}
// Amazing
template<class S>
void put_att( int ncid, int varid, const std::pair<S, nc_att_t>& att)
{
    S name = att.first;
    const nc_att_t& v = att.second;
    std::visit( [ncid, varid, name]( auto&& arg) { put_att( ncid, varid,
                std::make_pair( name, arg)); }, v);
}


template<class Iterable> // *it must be usable in put_att
void put_atts( int ncid, int varid, const Iterable& atts)
{
    for( const auto& it : atts)
        put_att( ncid, varid, it);
}

/////////////////////////////GETTERS////////////////////////////

template<class T>
std::vector<T> get_att_v( int ncid, int varid, std::string att)
{
    auto name = att.c_str();
    size_t size;
    NC_Error_Handle err;
    err = nc_inq_attlen( ncid, varid, name, &size);
    nc_type xtype;
    err = nc_inq_atttype( ncid, varid, name, &xtype);
    if constexpr( std::is_same_v<T, std::string>)
    {
        if( xtype == NC_STRING)
        {
            // size is the number of strings
            std::vector<char*> strs( size);
            err = nc_get_att_string( ncid, varid, name, strs.data());
            std::vector<T> data( strs.begin(), strs.end());
            nc_free_string( size, strs.data());
            return data;
        }
        // size is the number of characters (NetCDF returns NC_ECHAR if the
        // attribute is not text)
        std::string str( size, 'x');
        err = nc_get_att_text( ncid, varid, name, str.data());
        return std::vector<T>( 1, str);
    }
    else if constexpr( std::is_same_v<T, bool>) // read as NC_BYTE, 0 is false
    {
        // std::vector<bool> is not necessarily contiguous
        std::vector<signed char> tmp( size);
        err = nc_get_att_schar( ncid, varid, name, tmp.data());
        return std::vector<T>( tmp.begin(), tmp.end());
    }
    else if constexpr( is_nc_att_type_v<T>)
    {
        std::vector<T> data( size);
        err = get_att_T( ncid, varid, name, data.data());
        return data;
    }
    else // user defined type (e.g. compound), read without conversion
    {
        size_t type_size;
        err = nc_inq_type( ncid, xtype, NULL, &type_size);
        if( type_size != sizeof(T))
            throw std::runtime_error( "Size of given type does not match "
                "NetCDF type of attribute "+att);
        std::vector<T> data( size);
        err = nc_get_att( ncid, varid, name, data.data());
        return data;
    }
}

inline dg::file::nc_att_t get_att_t( int ncid, int varid,
        std::string att_name)
{
    auto name = att_name.c_str();
    size_t size;
    NC_Error_Handle err;
    err = nc_inq_attlen( ncid, varid, name, &size);
    nc_type xtype;
    err = nc_inq_atttype( ncid, varid, name, &xtype);
    // Signed integer types are read as int, unsigned integer types as unsigned
    // (NetCDF returns NC_ERANGE if 64 bit values do not fit)
    bool is_int = ( xtype == NC_SHORT or xtype == NC_INT or xtype == NC_INT64);
    bool is_uint = ( xtype == NC_UBYTE or xtype == NC_USHORT or xtype ==
        NC_UINT or xtype == NC_UINT64);
    if ( is_int and size == 1)
        return get_att_v<int>( ncid, varid, att_name)[0];
    else if ( is_int and size != 1)
        return get_att_v<int>( ncid, varid, att_name);
    else if ( is_uint and size == 1)
        return get_att_v<unsigned>( ncid, varid, att_name)[0];
    else if ( is_uint and size != 1)
        return get_att_v<unsigned>( ncid, varid, att_name);
    else if ( xtype == NC_FLOAT and size == 1)
        return get_att_v<float>( ncid, varid, att_name)[0];
    else if ( xtype == NC_FLOAT and size != 1)
        return get_att_v<float>( ncid, varid, att_name);
    else if ( xtype == NC_DOUBLE and size == 1)
        return get_att_v<double>( ncid, varid, att_name)[0];
    else if ( xtype == NC_DOUBLE and size != 1)
        return get_att_v<double>( ncid, varid, att_name);
    else if ( xtype == NC_BYTE and size == 1)
    {
        // BugFix: explicitly convert to bool
        bool value = get_att_v<bool>( ncid, varid, att_name)[0];
        return value;
    }
    else if ( xtype == NC_BYTE and size != 1)
        return get_att_v<bool>( ncid, varid, att_name);
    else if ( xtype == NC_CHAR or (xtype == NC_STRING and size == 1))
        return get_att_v<std::string>( ncid, varid, att_name)[0];
    else
        throw std::runtime_error( "Cannot convert attribute "+att_name+" to nc_att_t");
}

//namespace detail
//{
// utility overloads to be able to implement get_att and get_atts
template<class T>
void get_att_h( int ncid, int varid, std::string att_name, T& att)
{
    auto vec = get_att_v<T>( ncid, varid, att_name);
    if( vec.empty())
        throw std::runtime_error( "Attribute "+att_name+" is empty");
    att = vec[0];
}
template<class T>
void get_att_h( int ncid, int varid, std::string att_name, std::vector<T>& att)
{
    att = get_att_v<T>( ncid, varid, att_name);
}
inline void get_att_h( int ncid, int varid, std::string att_name, dg::file::nc_att_t& att)
{
    att = get_att_t( ncid, varid, att_name);
}
//}
template<class T> // T can be nc_att_t
T get_att_as( int ncid, int varid, std::string att_name)
{
    T att;
    detail::get_att_h( ncid, varid, att_name, att);
    return att;
}
template<class T> // T can be nc_att_t
std::map<std::string, T> get_atts_as( int ncid, int varid)
{
    NC_Error_Handle err;
    int number;
    if( varid == NC_GLOBAL)
        err = nc_inq_natts( ncid, &number);
    else
        err = nc_inq_varnatts( ncid, varid, &number);
    std::map < std::string, T> map;
    for( int i=0; i<number; i++)
    {
        char att_name[NC_MAX_NAME+1]; // 256 + null
        err = nc_inq_attname( ncid, varid, i, att_name);
        detail::get_att_h( ncid, varid, att_name, map[att_name]);
    }
    return map;
}
} // namespace detail
///@endcond
}// namespace file
}// namespace dg
