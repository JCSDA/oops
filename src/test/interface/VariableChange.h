/*
 * (C) Copyright 2018  UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_VARIABLECHANGE_H_
#define TEST_INTERFACE_VARIABLECHANGE_H_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/VariableChangeBase.h"
#include "oops/generic/instantiateVariableChangeFactories.h"
#include "oops/interface/State.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> class VariableChangeFixture : private boost::noncopyable {
  typedef oops::State<MODEL>                     State_;

 public:
  static std::vector<eckit::LocalConfiguration> & confs() {return getInstance().confs_;}
  static const State_      & xx()               {return *getInstance().xx_;}

 private:
  static VariableChangeFixture<MODEL>& getInstance() {
    static VariableChangeFixture<MODEL> theVariableChangeFixture;
    return theVariableChangeFixture;
  }

  VariableChangeFixture<MODEL>() {
    oops::instantiateVariableChangeFactories<MODEL>();

    TestEnvironment::config().get("VariableChangeTests", confs_);
  }

  ~VariableChangeFixture<MODEL>() {}

  std::vector<eckit::LocalConfiguration>             confs_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testVariableChangeInverse() {
}

// -----------------------------------------------------------------------------

template <typename MODEL> class VariableChange : public oops::Test {
 public:
  VariableChange() {}
  virtual ~VariableChange() {}
 private:
  std::string testid() const {return "test::VariableChange<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/VariableChange");

    ts->add(BOOST_TEST_CASE(&testVariableChangeInverse<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_VARIABLECHANGE_H_
