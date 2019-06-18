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
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class ObsAuxControlFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & config() {return *getInstance().conf_;}

 private:
  static ObsAuxControlFixture<MODEL>& getInstance() {
    static ObsAuxControlFixture<MODEL> theObsAuxControlFixture;
    return theObsAuxControlFixture;
  }

  ObsAuxControlFixture() {
    std::vector<eckit::LocalConfiguration> osconf;
    TestEnvironment::config().get("Observations.ObsTypes", osconf);
    conf_.reset(new eckit::LocalConfiguration(osconf[0]));
  }

  ~ObsAuxControlFixture() {}

  std::unique_ptr<const eckit::LocalConfiguration>  conf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsAuxControlFixture<MODEL>   Test_;
  typedef oops::ObsAuxControl<MODEL>    ObsAux_;

  std::unique_ptr<ObsAux_> bias(new ObsAux_(Test_::config()));
  EXPECT(bias.get());

  bias.reset();
  EXPECT(!bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ObsAuxControlFixture<MODEL>   Test_;
  typedef oops::ObsAuxControl<MODEL>    ObsAux_;

  std::unique_ptr<ObsAux_> bias(new ObsAux_(Test_::config()));

  std::unique_ptr<ObsAux_> other(new ObsAux_(*bias));
  EXPECT(other.get());

  other.reset();
  EXPECT(!other.get());

  EXPECT(bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsAuxControl : public oops::Test {
 public:
  ObsAuxControl() {}
  virtual ~ObsAuxControl() {}
 private:
  std::string testid() const {return "test::ObsAuxControl<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsAuxControl/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsAuxControl/testCopyConstructor")
      { testCopyConstructor<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXCONTROL_H_
