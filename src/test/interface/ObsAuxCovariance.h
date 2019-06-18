/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSAUXCOVARIANCE_H_
#define TEST_INTERFACE_OBSAUXCOVARIANCE_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/interface/ObsAuxCovariance.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class ObsAuxCovarianceFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & config() {return *getInstance().conf_;}

 private:
  static ObsAuxCovarianceFixture<MODEL>& getInstance() {
    static ObsAuxCovarianceFixture<MODEL> theObsAuxCovarianceFixture;
    return theObsAuxCovarianceFixture;
  }

  ObsAuxCovarianceFixture() {
    std::vector<eckit::LocalConfiguration> osconf;
    TestEnvironment::config().get("Observations.ObsTypes", osconf);
    conf_.reset(new eckit::LocalConfiguration(osconf[0]));
  }

  ~ObsAuxCovarianceFixture() {}

  std::unique_ptr<const eckit::LocalConfiguration>  conf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsAuxCovarianceFixture<MODEL>   Test_;
  typedef oops::ObsAuxCovariance<MODEL>    Covariance_;

  std::unique_ptr<Covariance_> bias(new Covariance_(Test_::config()));
  EXPECT(bias.get());

  bias.reset();
  EXPECT(!bias.get());
}

// -----------------------------------------------------------------------------
//  void linearize(const ObsAuxControl_ &, const Geometry_ &);
//  void multiply(const ObsAuxIncrement_ &, ObsAuxIncrement_ &) const;
//  void invMult(const ObsAuxIncrement_ &, ObsAuxIncrement_ &) const;
//  void randomize(ObsAuxIncrement_ &) const;
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsAuxCovariance : public oops::Test {
 public:
  ObsAuxCovariance() {}
  virtual ~ObsAuxCovariance() {}
 private:
  std::string testid() const {return "test::ObsAuxCovariance<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsAuxCovariance/testConstructor")
      { testConstructor<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXCOVARIANCE_H_
