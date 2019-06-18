/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_BASE_OBSERRORCOVARIANCE_H_
#define TEST_BASE_OBSERRORCOVARIANCE_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0


#include "eckit/testing/Test.h"
#include "oops/base/ObsErrorBase.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsTestsFixture<MODEL>  Test_;
  typedef oops::ObsErrorBase<MODEL>  Covar_;
  typedef oops::ObsVector<MODEL>     ObsVector_;

  oops::instantiateObsErrorFactory<MODEL>();

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ obserr(Test_::obspace()[jj], "ObsError");
    obserr.save("EffectiveError");

    const eckit::LocalConfiguration rconf(conf[jj], "Covariance");
    std::unique_ptr<Covar_> R(
      oops::ObsErrorFactory<MODEL>::create(rconf, Test_::obspace()[jj]));
    EXPECT(R.get());

    R.reset();
    EXPECT(!R.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsErrorCovariance : public oops::Test {
 public:
  ObsErrorCovariance() {}
  virtual ~ObsErrorCovariance() {}
 private:
  std::string testid() const {return "test::ObsErrorCovariance<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsErrorCovariance/testConstructor")
      { testConstructor<MODEL>(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_OBSERRORCOVARIANCE_H_
