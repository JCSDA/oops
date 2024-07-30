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
/// \brief test constructor and print method
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsAuxControl<OBS>    ObsAux_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
    std::unique_ptr<ObsAux_> bias(new ObsAux_(Test_::obspace()[jj], biasconf));
    EXPECT(bias.get());
    oops::Log::info() << "Testing ObsAuxControl: " << *bias << std::endl;

    // Not all configurations for interface tests specify "obs bias"; need to check
    // whether "obs bias" section is available
    if (Test_::config(jj).has("obs bias")) {
      const double reference = Test_::config(jj).getDouble("obs bias test.norm");
      const double tolerance = Test_::config(jj).getDouble("obs bias test.relative tolerance");
      EXPECT(oops::is_close_relative(bias->norm(), reference, tolerance));
    }

    bias.reset();
    EXPECT(!bias.get());
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testCopyConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsAuxControl<OBS>    ObsAux_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
    std::unique_ptr<ObsAux_> bias(new ObsAux_(Test_::obspace()[jj], biasconf));

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
  typedef ObsTestsFixture<OBS> Test_;
 public:
  using oops::Test::Test;
  virtual ~ObsAuxControl() {}
 private:
  std::string testid() const override {return "test::ObsAuxControl<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsAuxControl/testConstructor")
      { testConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsAuxControl/testCopyConstructor")
      { testCopyConstructor<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXCONTROL_H_
