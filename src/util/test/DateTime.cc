/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include "util/DateTime.h"
#include "util/Duration.h"

#include <string>

#include <boost/test/unit_test.hpp>

namespace {

BOOST_AUTO_TEST_SUITE(test_DateTime)

// -----------------------------------------------------------------------------

  BOOST_AUTO_TEST_CASE(test_construct_from_string) {
    std::string s1("2011-09-16T13:55:20Z");
    util::DateTime d1(s1);
    std::string s2(d1.toString());
    BOOST_CHECK_EQUAL(s2, s1);
  }

// -----------------------------------------------------------------------------

//  This method is now private. Can we test private methods?

//  BOOST_AUTO_TEST_CASE(test_set_from_string) {
//    util::DateTime d1("2012-12-05T15:42:55Z");
//    std::string s1("2011-09-16T13:55:20Z");
//    d1.set(s1);
//    std::string s2(d1.toString());
//    BOOST_CHECK_EQUAL(s2, s1);
//  }

// -----------------------------------------------------------------------------

  BOOST_AUTO_TEST_CASE(test_copy_constructor) {
    std::string s1("2011-09-16T13:55:20Z");
    util::DateTime d1(s1);
    util::DateTime d2(d1);
    std::string s2(d2.toString());
    BOOST_CHECK_EQUAL(s2, s1);
  }

// -----------------------------------------------------------------------------

  BOOST_AUTO_TEST_CASE(test_operators) {
     util::DateTime d1(2011, 9, 16, 13, 55, 20);
     util::DateTime d2(2011, 9, 16, 13, 55, 21);
     util::DateTime d3(2011, 9, 16, 13, 55, 21);
     BOOST_CHECK(d2 == d3);
     BOOST_CHECK(d1 != d2);
     BOOST_CHECK(d1 <  d2);
     BOOST_CHECK(d3 <= d2);
     BOOST_CHECK(d2 >  d1);
     BOOST_CHECK(d3 >= d2);

     util::Duration dur(d2-d1);  // one second

     BOOST_CHECK_EQUAL((d1+dur), d2);
     BOOST_CHECK_EQUAL((d2-dur), d1);

     util::DateTime d4(d1);
     d4 += dur;
     BOOST_CHECK_EQUAL(d4, d2);
     d4 -= dur;
     BOOST_CHECK_EQUAL(d4, d1);
  }

// -----------------------------------------------------------------------------

  BOOST_AUTO_TEST_CASE(test_stream_write) {
    std::string s1("2011-09-16T13:55:20Z");
    util::DateTime d1(s1);
    std::stringstream strm;
    strm << d1;
    std::string s2;
    strm >> s2;
    BOOST_CHECK_EQUAL(s2, s1);
  }

// -----------------------------------------------------------------------------

  BOOST_AUTO_TEST_CASE(test_stream_read) {
    std::string s1("2011-09-16T13:55:20Z");
    std::stringstream strm;
    strm << s1;

    util::DateTime d1("2012-12-05T15:42:55Z");
    strm >> d1;
    BOOST_CHECK_EQUAL(s1, d1.toString());
  }

// -----------------------------------------------------------------------------

  BOOST_AUTO_TEST_CASE(test_special_dates) {
    //  Y2K tests
    util::Duration oneDay("P1D");  // one day
    util::DateTime d1(2000, 2, 28, 0, 0, 0);
    BOOST_CHECK_EQUAL((d1+oneDay).toString(), "2000-02-29T00:00:00Z");

    util::DateTime d2(2000, 2, 29, 0, 0, 0);
    BOOST_CHECK_EQUAL((d2-oneDay).toString(), "2000-02-28T00:00:00Z");

    util::DateTime d3(2000, 3, 1, 0, 0, 0);
    BOOST_CHECK_EQUAL((d3-oneDay).toString(), "2000-02-29T00:00:00Z");

    // Far in the past/future
    util::Duration longDur("P345126D");
    util::DateTime d4(2011, 9, 16, 9, 0, 0);
    BOOST_CHECK_EQUAL((d4 - longDur).toString(), "1066-10-14T09:00:00Z");
    BOOST_CHECK_EQUAL((d4 + longDur).toString(), "2956-08-18T09:00:00Z");
  }

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
