/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/testing/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace {


// -----------------------------------------------------------------------------

  CASE("test_construct_from_int") {
    util::Duration dur(86400);
    EXPECT(dur.toSeconds() == 86400);
  }

// -----------------------------------------------------------------------------

  CASE("test_construct_from_string") {
    util::Duration dur("-P1DT2H5M39S");
    EXPECT(dur.toSeconds() == -93939);
  }

// -----------------------------------------------------------------------------

//  Private method: how do we test?

//  CASE("test_set_from_string") {
//    util::Duration dur(0);
//    dur.set("P1DT2H5M39S");
//    EXPECT(dur.toSeconds() == 93939);
//
//    dur.set("-P1D");
//    EXPECT(dur.toSeconds() == -86400);
//    dur.set("PT3H");
//    EXPECT(dur.toSeconds() == 10800);
//    dur.set("-PT5M");
//    EXPECT(dur.toSeconds() == -300);
//    dur.set("PT59S");
//    EXPECT(dur.toSeconds() == 59);
//  }

// -----------------------------------------------------------------------------

  CASE("test_convert_to_string") {
    util::Duration dur(93939);
    std::string s("P1DT2H5M39S");
    EXPECT(dur.toString() == s);
  }

// -----------------------------------------------------------------------------

  CASE("test_operators") {
    util::Duration dur0(0);
    util::Duration dur1(1);
    util::Duration dur1a(1);

    EXPECT(dur1 == dur1a);
    EXPECT(dur0 != dur1);
    EXPECT(dur0 <  dur1);
    EXPECT(dur1 <= dur1a);
    EXPECT(dur1 >  dur0);
    EXPECT(dur1 >= dur1a);

    util::Duration dur1b(-1);
    dur1b.negate();
    EXPECT(dur1b == dur1);

    util::Duration dur1c(-1);
    EXPECT(-dur1c == dur1);

    util::Duration dur(86400);
    dur -= dur1;
    EXPECT(dur.toSeconds() == 86399);
    dur += dur1;
    EXPECT(dur.toSeconds() == 86400);

    util::Duration dayAndABit(86499);
    util::Duration hour(3600);
    EXPECT(dayAndABit%hour == 99);
  }

// -----------------------------------------------------------------------------

  CASE("test_long_durations") {
    util::Duration dur1("P4000000DT1S");  // about 10000 years
    util::Duration dur2("-P4000000D");
    dur1 += dur2;
    EXPECT(dur1.toSeconds() == 1);
  }

// -----------------------------------------------------------------------------

  CASE("test_stream_write") {
    std::string s1("P1DT2H5M39S");
    util::Duration d1(s1);
    std::stringstream strm;
    strm << d1;
    std::string s2;
    strm >> s2;
    EXPECT(s2 == s1);
  }

// -----------------------------------------------------------------------------

  CASE("test_stream_read") {
    std::string s1("P1DT2H5M39S");
    std::stringstream strm;
    strm << s1;
    util::Duration d1;
    strm >> d1;
    EXPECT(s1 == d1.toString());
  }

// -----------------------------------------------------------------------------

}  // anonymous namespace

int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
