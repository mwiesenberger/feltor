#include <iostream>
#include <cassert>
#include <filesystem>

#include "catch2/catch_all.hpp"
#include "json_wrapper.h"

// TODO Needs to be more precise
#ifdef DG_USE_JSONHPP
TEST_CASE("Json HPP wrapper")
#else
TEST_CASE("Json CPP wrapper")
#endif
{
    std::ofstream file( "json_wrapper.json");
    file << R"unique({
        "hello": "world",
        "array": [ 42, 73],
        "nested":
        {
            "hi" : 38.0,
            "bool" : true,
            "another_nest":
            {
                "number" : 7
            }
        }
    })unique";
    file.close();
    REQUIRE( std::filesystem::exists( "json_wrapper.json"));
    auto js = dg::file::file2Json( "json_wrapper.json");
    std::string hello0, hello1;
    int idx0, idx1, idx2, idx3;
    unsigned uint0, uint1, uint2;
    double hi;
    bool test;
    SECTION( "Silent")
    {
        dg::file::WrappedJsonValue ws( js, dg::file::error::is_silent);
        std::stringstream ss;
        // catch stream buffer
        auto cerr_buf = std::cerr.rdbuf( ss.rdbuf() );
        hello0 = ws[ "hell"].asString( "default" );
        hello1 = ws.get( "hell", "default").asString( );
        idx0 = ws[ "array"][2].asInt();
        idx1 = ws[ "array"].get(2,0).asInt();
        idx2 = ws["hello"].asInt();
        idx3 = ws[ "arr"][1].asInt( );
        hi = ws[ "neted"]["hi"].asDouble( );
        test = ws[ "nested"]["bol"].asBool( true);
        uint0 = ws["nested"]["another_net"]["number"].asUInt();
        uint1 = ws["hello"]["another_net"]["number"].asUInt();
        uint2 = ws[0][2]["number"].asUInt();
        // reset buffer
        std::cerr.rdbuf( cerr_buf );
        CHECK( ss.str().empty());
        //Avoid compiler warnings
        ss << idx0<<idx1<<idx2<<idx3<<hi<<test<<uint0<<uint1<<uint2;
    }
    SECTION( "Warning")
    {
        dg::file::WrappedJsonValue ws( js, dg::file::error::is_warning);
        std::stringstream ss;
        // catch stream buffer
        auto cerr_buf = std::cerr.rdbuf( ss.rdbuf() );
        hello0 = ws[ "hell"].asString( "default" );
        CHECK( not ss.str().empty());
        ss.str( "");
        hello1 = ws.get( "hell", "default").asString( );
        CHECK( not ss.str().empty());
        ss.str( "");
        idx0 = ws[ "array"][2].asInt();
        CHECK( not ss.str().empty());
        ss.str( "");
        idx1 = ws[ "array"].get(2,0).asInt();
        CHECK( not ss.str().empty());
        ss.str( "");
        idx2 = ws["hello"].asInt();
        CHECK( not ss.str().empty());
        ss.str( "");
        idx3 = ws[ "arr"][1].asInt( );
        CHECK( not ss.str().empty());
        ss.str( "");
        hi = ws[ "neted"]["hi"].asDouble( );
        CHECK( not ss.str().empty());
        ss.str( "");
        test = ws[ "nested"]["bol"].asBool( true);
        CHECK( not ss.str().empty());
        ss.str( "");
        uint0 = ws["nested"]["another_net"]["number"].asUInt();
        CHECK( not ss.str().empty());
        ss.str( "");
        uint1 = ws["hello"]["another_net"]["number"].asUInt();
        CHECK( not ss.str().empty());
        ss.str( "");
        uint2 = ws[0][2]["number"].asUInt();
        CHECK( not ss.str().empty());
        ss.str( "");
        // reset buffer
        std::cerr.rdbuf( cerr_buf );
        //Avoid compiler warnings
        ss << idx0<<idx1<<idx2<<idx3<<hi<<test<<uint0<<uint1<<uint2;
    }
    SECTION( "Error throws")
    {
        // We test that it throws and that the message contains certain key words
        using Catch::Matchers::ContainsSubstring;
        dg::file::WrappedJsonValue ws( js, dg::file::error::is_throw);
        CHECK_THROWS_WITH( ws[ "hell"].asString( "default" ),
            ContainsSubstring("\"hell\""));
        CHECK_THROWS_WITH( ws.get( "hell", "default").asString( ),
            ContainsSubstring("\"hell\""));
        CHECK_THROWS_WITH( ws[ "array"][2].asInt(),
            ContainsSubstring("\"array\"") && ContainsSubstring("2"));
        CHECK_THROWS_WITH( ws[ "array"].get(2,0).asInt(),
            ContainsSubstring("\"array\"") && ContainsSubstring("2"));
        CHECK_THROWS_WITH( ws["hello"].asInt(),
            ContainsSubstring("\"hello\"") && ContainsSubstring("\"world\"")
            && ContainsSubstring("Int"));
        CHECK_THROWS_WITH( ws[ "arr"][1].asInt( ),
            ContainsSubstring("\"arr\""));
        CHECK_THROWS_WITH( ws[ "neted"]["hi"].asDouble( ),
            ContainsSubstring("\"neted\""));
        CHECK_THROWS_WITH( ws[ "nested"]["bol"].asBool( true),
            ContainsSubstring("\"nested\"") && ContainsSubstring("\"bol\""));
        CHECK_THROWS_WITH( ws["nested"]["another_net"]["number"].asUInt(),
            ContainsSubstring("\"nested\"") && ContainsSubstring("\"another_net\""));
        CHECK_THROWS_WITH( ws["hello"]["another_net"]["number"].asUInt(),
            ContainsSubstring("\"hello\"") && ContainsSubstring("\"another_net\""));
        CHECK_THROWS_WITH( ws[0][2]["number"].asUInt(),
            ContainsSubstring("0"));
    }
    SECTION("Access")
    {
        dg::file::WrappedJsonValue ws( js);
        std::string hello = ws["hello"].asString();
        CHECK( hello == "world");
        idx0 = ws[ "array"][0].asInt();
        CHECK( ws["array"].size() == 2);
        CHECK( idx0 == 42);
        idx0 = ws[ "array"].get(0,0).asInt();
        CHECK( idx0 == 42);
        idx1 = ws["array"][1].asInt();
        CHECK( idx1 == 73);
        hi = ws[ "nested"][ "hi"].asDouble();
        CHECK( hi == 38);
        test = ws[ "nested"]["bool"].asBool(true);
        CHECK( test == true);
        unsigned uint = ws["nested"]["another_nest"]["number"].asUInt();
        CHECK( uint == 7);
        INFO( "Test access string");
        std::string acc = ws["nested"]["another_nest"]["number"].access_string();
        INFO( acc);
        CHECK( acc.find("\"nested\"") != std::string::npos);
        CHECK( acc.find("\"another_nest\"") != std::string::npos);
        CHECK( acc.find("\"number\"") != std::string::npos);
    }
    SECTION("Existence")
    {
        dg::file::WrappedJsonValue ws( js);
        REQUIRE( ws.isMember( "hello"));
        CHECK( ws["hello"].isString());
        CHECK( not ws["array"].isString());
        CHECK( not ws["nested"].isString());
        CHECK( not ws["nested"]["hi"].isString());
        CHECK( not ws["nested"]["bool"].isString());
    }

    SECTION( "Test Vector 2 Json")
    {
        dg::file::WrappedJsonValue ws( dg::file::vec2json( {42, 73}));
        CHECK( ws.get( 0,0).asInt() == 42);
        CHECK( ws.get( 1,0).asInt() == 73);
    }
    SECTION( "Size of values")
    {
        dg::file::WrappedJsonValue ws( js);
        CHECK( ws["array"].size() == 2);
        CHECK( ws["nested"].size() == 3);
        CHECK( ws["hello"].size() == 0);
        CHECK( ws["nested"]["hi"].size() == 0);
    }
    SECTION( "Integers are checked strictly")
    {
        dg::file::WrappedJsonValue ws( dg::file::string2Json( R"({
            "neg" : -1, "frac" : 2.5, "integral" : 3.0, "big" : 1e20,
            "max" : 4294967295})"));
        CHECK( ws["neg"].asInt() == -1);
        CHECK_THROWS_AS( ws["neg"].asUInt(), std::runtime_error);
        CHECK_THROWS_AS( ws["frac"].asInt(), std::runtime_error);
        CHECK_THROWS_AS( ws["frac"].asUInt(), std::runtime_error);
        CHECK( ws["frac"].asDouble() == 2.5);
        CHECK( ws["integral"].asInt() == 3);
        CHECK( ws["integral"].asUInt() == 3);
        CHECK_THROWS_AS( ws["big"].asInt(), std::runtime_error);
        CHECK_THROWS_AS( ws["big"].asUInt(), std::runtime_error);
        CHECK( ws["max"].asUInt() == 4294967295u);
        CHECK_THROWS_AS( ws["max"].asInt(), std::runtime_error);
        ws.set_mode( dg::file::error::is_silent);
        CHECK( ws["neg"].asUInt( 7) == 7);
    }
    SECTION( "Comment handling")
    {
        // jsoncpp does not forbid comments directly after "{" or before "}"
        std::string commented = R"({
            "a" : // This is a comment
            1
        })";
        for( auto comm : {dg::file::comments::are_kept,
                          dg::file::comments::are_discarded})
        {
            dg::file::WrappedJsonValue ws( dg::file::string2Json( commented,
                comm));
            CHECK( ws["a"].asInt() == 1);
        }
        CHECK_THROWS_AS( dg::file::string2Json( commented,
            dg::file::comments::are_forbidden), std::runtime_error);
    }
    SECTION( "Parse error handling")
    {
        std::string invalid = R"({ "a" : })";
        CHECK_THROWS_AS( dg::file::string2Json( invalid), std::runtime_error);
        CHECK_THROWS_AS( dg::file::file2Json( "does-not-exist.json"),
            std::runtime_error);

        std::stringstream ss;
        auto cerr_buf = std::cerr.rdbuf( ss.rdbuf() );
        auto js_silent = dg::file::string2Json( invalid,
            dg::file::comments::are_discarded, dg::file::error::is_silent);
        auto file_silent = dg::file::file2Json( "does-not-exist.json",
            dg::file::comments::are_discarded, dg::file::error::is_silent);
        bool silent = ss.str().empty();
        auto js_warning = dg::file::string2Json( invalid,
            dg::file::comments::are_discarded, dg::file::error::is_warning);
        bool warning = not ss.str().empty();
        std::cerr.rdbuf( cerr_buf );
        CHECK( silent);
        CHECK( warning);
        CHECK( dg::file::WrappedJsonValue( js_silent).size() == 0);
        CHECK( dg::file::WrappedJsonValue( js_warning).size() == 0);
        CHECK( dg::file::WrappedJsonValue( file_silent).size() == 0);
    }
    std::filesystem::remove( "json_wrapper.json");

}
