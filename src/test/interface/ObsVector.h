/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSVECTOR_H_
#define TEST_INTERFACE_OBSVECTOR_H_

#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/ObsVector.h"
#include "test/TestEnvironment.h"
#include "test/interface/ObsTestsFixture.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsTestsFixture<MODEL>  Test_;
  typedef oops::ObsVector<MODEL>  ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<ObsVector_> ov(new ObsVector_(Test_::obspace()[jj]));
    BOOST_CHECK(ov.get());

    ov.reset();
    BOOST_CHECK(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ObsTestsFixture<MODEL>  Test_;
  typedef oops::ObsVector<MODEL>  ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<ObsVector_> ov(new ObsVector_(Test_::obspace()[jj]));

    boost::scoped_ptr<ObsVector_> other(new ObsVector_(*ov));
    BOOST_CHECK(other.get());

    other.reset();
    BOOST_CHECK(!other.get());

    BOOST_CHECK(ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ObsVector : public oops::Test {
 public:
  ObsVector() {}
  virtual ~ObsVector() {}
 private:
  std::string testid() const {return "test::ObsVector<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObsVector");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testCopyConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSVECTOR_H_
