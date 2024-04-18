/*
 * (C) Copyright 2009-2020 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/DateTime.h"

#include <sstream>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "eckit/testing/Test.h"

#include "oops/util/Duration.h"

namespace {


// -----------------------------------------------------------------------------

  CASE("test_construct_from_string") {
    std::string s1("2011-09-16T13:55:20Z");
    util::DateTime d1(s1);
    std::string s2(d1.toString());
    EXPECT(s2 == s1);
  }

// -----------------------------------------------------------------------------

  CASE("test_construct_from_ymdhms") {
    std::string s1("2011-09-16T13:55:20Z");
    util::DateTime d1(2011, 9, 16, 13, 55, 20);
    std::string s2(d1.toString());
    EXPECT(s2 == s1);
  }

// -----------------------------------------------------------------------------

  CASE("test_construct_from_YYYYMMDD_hhmmss") {
    std::string s1("2011-09-16T13:55:20Z");
    util::DateTime d1(20110916, 135520);
    std::string s2(d1.toString());
    EXPECT(s2 == s1);
  }

// -----------------------------------------------------------------------------

//  This method is now private. Can we test private methods?

//  CASE("test_set_from_string") {
//    util::DateTime d1("2012-12-05T15:42:55Z");
//    std::string s1("2011-09-16T13:55:20Z");
//    d1.set(s1);
//    std::string s2(d1.toString());
//    EXPECT(s2 == s1);
//  }

// -----------------------------------------------------------------------------

  CASE("test_copy_constructor") {
    std::string s1("2011-09-16T13:55:20Z");
    util::DateTime d1(s1);
    util::DateTime d2(d1);
    std::string s2(d2.toString());
    EXPECT(s2 == s1);
  }

// -----------------------------------------------------------------------------

  CASE("test_uninitialized_date") {
    util::DateTime d1;
    EXPECT_THROWS_AS(d1.toString(), eckit::BadValue);
  }

// -----------------------------------------------------------------------------

  CASE("test_to_YYYYMMDD_hhmmss") {
    int date1 = 20200108;
    int time1 = 32745;
    util::DateTime d1(date1, time1);
    int date2;
    int time2;
    d1.toYYYYMMDDhhmmss(date2, time2);
    EXPECT(date2 == date1);
    EXPECT(time2 == time1);
  }

// -----------------------------------------------------------------------------

  CASE("test_seconds_since_jan1") {
    const int date = 20200518;
    const int time = 32745;
    const util::DateTime datetime(date, time);
    const int seconds_since_jan1 = datetime.secondsSinceJan1();
    // Below, the number 138 comes from the number of elapsed days before May 18th,
    // which was the 139th day of 2020
    const int expected = 45 + 27 * 60 + 3 * 3600 + 138 * 86400;
    EXPECT(seconds_since_jan1 == expected);
  }

// -----------------------------------------------------------------------------

  CASE("test_operators") {
     util::DateTime d1(2011, 9, 16, 13, 55, 20);
     util::DateTime d2(2011, 9, 16, 13, 55, 21);
     util::DateTime d3(2011, 9, 16, 13, 55, 21);
     EXPECT(d2 == d3);
     EXPECT(d1 != d2);
     EXPECT(d1 <  d2);
     EXPECT(d3 <= d2);
     EXPECT(d2 >  d1);
     EXPECT(d3 >= d2);

     util::Duration dur(d2-d1);  // one second

     EXPECT((d1+dur) == d2);
     EXPECT((d2-dur) == d1);

     util::DateTime d4(d1);
     d4 += dur;
     EXPECT(d4 == d2);
     d4 -= dur;
     EXPECT(d4 == d1);
  }

// -----------------------------------------------------------------------------

  CASE("test_stream_write") {
    std::string s1("2011-09-16T13:55:20Z");
    util::DateTime d1(s1);
    std::stringstream strm;
    strm << d1;
    std::string s2;
    strm >> s2;
    EXPECT(s2 == s1);
  }

// -----------------------------------------------------------------------------

  CASE("test_stream_read") {
    std::string s1("2011-09-16T13:55:20Z");
    std::stringstream strm;
    strm << s1;

    util::DateTime d1("2012-12-05T15:42:55Z");
    strm >> d1;
    EXPECT(s1 == d1.toString());
  }

// -----------------------------------------------------------------------------

  CASE("test_special_dates") {
    //  Y2K tests
    util::Duration oneDay("P1D");  // one day
    util::DateTime d1(2000, 2, 28, 0, 0, 0);
    EXPECT((d1+oneDay).toString() == "2000-02-29T00:00:00Z");

    util::DateTime d2(2000, 2, 29, 0, 0, 0);
    EXPECT((d2-oneDay).toString() == "2000-02-28T00:00:00Z");

    util::DateTime d3(2000, 3, 1, 0, 0, 0);
    EXPECT((d3-oneDay).toString() == "2000-02-29T00:00:00Z");

    // Far in the past/future
    util::Duration longDur("P345126D");
    util::DateTime d4(2011, 9, 16, 9, 0, 0);
    EXPECT((d4 - longDur).toString() == "1066-10-14T09:00:00Z");
    EXPECT((d4 + longDur).toString() == "2956-08-18T09:00:00Z");
  }

// -----------------------------------------------------------------------------

}  // anonymous namespace

int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
