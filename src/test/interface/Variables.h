/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_VARIABLES_H_
#define TEST_INTERFACE_VARIABLES_H_

#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/Variables.h"
#include "test/TestEnvironment.h"
#include "eckit/config/Configuration.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class VariablesFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & getConfig() {return *getInstance().conf_;}

 private:
  static VariablesFixture<MODEL>& getInstance() {
    static VariablesFixture<MODEL> theVariablesFixture;
    return theVariablesFixture;
  }

  VariablesFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "Variables"));
  }

  ~VariablesFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
};
// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef oops::Variables<MODEL>        Variables_;

  boost::scoped_ptr<Variables_> vars(new Variables_(VariablesFixture<MODEL>::getConfig()));
  BOOST_CHECK(vars.get());

  vars.reset();
  BOOST_CHECK(!vars.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef oops::Variables<MODEL>        Variables_;
  boost::scoped_ptr<Variables_> vars(new Variables_(VariablesFixture<MODEL>::getConfig()));

  boost::scoped_ptr<Variables_> other(new Variables_(*vars));
  BOOST_CHECK(other.get());

  other.reset();
  BOOST_CHECK(!other.get());

  BOOST_CHECK(vars.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> class Variables : public oops::Test {
 public:
  Variables() {}
  virtual ~Variables() {}
 private:
  std::string testid() const {return "test::Variables<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/Variables");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testCopyConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_VARIABLES_H_
