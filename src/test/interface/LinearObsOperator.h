/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_LINEAROBSOPERATOR_H_
#define TEST_INTERFACE_LINEAROBSOPERATOR_H_

#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/LinearObsOperator.h"
#include "test/TestEnvironment.h"
#include "test/interface/ObsTestsFixture.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsTestsFixture<MODEL>  Test_;
  typedef oops::LinearObsOperator<MODEL>  LinearObsOperator_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<LinearObsOperator_> ov(new LinearObsOperator_(Test_::obspace()[jj]));
    BOOST_CHECK(ov.get());

    ov.reset();
    BOOST_CHECK(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class LinearObsOperator : public oops::Test {
 public:
  LinearObsOperator() {}
  virtual ~LinearObsOperator() {}
 private:
  std::string testid() const {return "test::LinearObsOperator<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/LinearObsOperator");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEAROBSOPERATOR_H_
