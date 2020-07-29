/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSAUXCONTROL_H_
#define TEST_INTERFACE_OBSAUXCONTROL_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsAuxControl<OBS>    ObsAux_;

  std::vector<eckit::LocalConfiguration> oconf;
  TestEnvironment::config().get("observations", oconf);
  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    std::unique_ptr<ObsAux_> bias(new ObsAux_(Test_::obspace()[jj], oconf[jj]));
    EXPECT(bias.get());

    bias.reset();
    EXPECT(!bias.get());
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testCopyConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsAuxControl<OBS>    ObsAux_;

  std::vector<eckit::LocalConfiguration> oconf;
  TestEnvironment::config().get("observations", oconf);
  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    std::unique_ptr<ObsAux_> bias(new ObsAux_(Test_::obspace()[jj], oconf[jj]));

    std::unique_ptr<ObsAux_> other(new ObsAux_(*bias));
    EXPECT(other.get());

    other.reset();
    EXPECT(!other.get());

    EXPECT(bias.get());
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsAuxControl : public oops::Test {
 public:
  ObsAuxControl() {}
  virtual ~ObsAuxControl() {}
 private:
  std::string testid() const {return "test::ObsAuxControl<" + OBS::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsAuxControl/testConstructor")
      { testConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsAuxControl/testCopyConstructor")
      { testCopyConstructor<OBS>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXCONTROL_H_
