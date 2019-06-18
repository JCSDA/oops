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

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"
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
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/VariableChange/testVariableChangeInverse")
      { testVariableChangeInverse<MODEL>(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_VARIABLECHANGE_H_
