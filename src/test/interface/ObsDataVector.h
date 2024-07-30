/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_OBSDATAVECTOR_H_
#define TEST_INTERFACE_OBSDATAVECTOR_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/runs/Test.h"
#include "oops/util/dot_product.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief tests constructor and print method
template <typename OBS> void testConstructors() {
  typedef ObsTestsFixture<OBS>             Test_;
  typedef oops::ObsDataVector<OBS, float>  ObsDataVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    oops::ObsVariables vars(Test_::obspace()[jj].obsvariables());
    std::unique_ptr<ObsDataVector_> odv1(new ObsDataVector_(Test_::obspace()[jj], vars));
    EXPECT(odv1.get());

    odv1->zero();
    oops::Log::info() << "Printing zero ObsDataVector: " << *odv1 << std::endl;

    std::unique_ptr<ObsDataVector_> odv2(new ObsDataVector_(*odv1));
    EXPECT(odv2.get());

    odv2.reset();
    EXPECT(!odv2.get());
    EXPECT(odv1.get());

    odv1.reset();
    EXPECT(!odv1.get());
  }
}

// -----------------------------------------------------------------------------
template <typename OBS> void testObsVector() {
  typedef ObsTestsFixture<OBS>             Test_;
  typedef oops::ObsDataVector<OBS, float>  ObsDataVector_;

  const double tolerance = 1.0e-10;
  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    oops::ObsVector<OBS> ov1(Test_::obspace()[jj]);
    ov1.random();

    ObsDataVector_ odv(ov1);
    oops::Log::info() << "Printing random ObsDataVector: " << odv << std::endl;

    oops::ObsVector<OBS> ov2(Test_::obspace()[jj]);
    ov2 = odv;

    ov1 -= ov2;
    oops::Log::info() << "ObsVector - ObsDataVector = " << ov1 << std::endl;
    const double diff = dot_product(ov1, ov1);
    oops::Log::info() << "ObsVector, ObsDataVector diff = " << diff << std::endl;
    EXPECT(diff < tolerance);
  }
}
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsDataVector : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  using oops::Test::Test;
  virtual ~ObsDataVector() {}

 private:
  std::string testid() const override {return "test::ObsDataVector<" + OBS::name() + ", float>";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsDataVector/testConstructors")
      { testConstructors<OBS>(); });
    ts.emplace_back(CASE("interface/ObsDataVector/testObsVector")
      { testObsVector<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSDATAVECTOR_H_
