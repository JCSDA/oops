/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_LOCALIZATION_H_
#define TEST_INTERFACE_LOCALIZATION_H_

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Localization.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> class LocalizationFixture : private boost::noncopyable {
  typedef oops::Localization<MODEL>                     Localization_;
  typedef oops::Geometry<MODEL>                         Geometry_;
  typedef oops::IncrementEnsemble<MODEL>                Ensemble_;

 public:
  static const Geometry_       & resol()        {return *getInstance().resol_;}
  static const oops::Variables & ctlvars()      {return *getInstance().ctlvars_;}
  static const util::DateTime  & time()         {return *getInstance().time_;}
  static const Localization_   & localization() {return *getInstance().local_;}

 private:
  static LocalizationFixture<MODEL>& getInstance() {
    static LocalizationFixture<MODEL> theLocalizationFixture;
    return theLocalizationFixture;
  }

  LocalizationFixture<MODEL>() {
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

    ctlvars_.reset(new oops::Variables(TestEnvironment::config(), "loc variables"));

    time_.reset(new util::DateTime(TestEnvironment::config().getString("test date")));

//  Setup the localization matrix
    const eckit::LocalConfiguration conf(TestEnvironment::config(), "localization");
    local_.reset(new Localization_(*resol_, conf));

    oops::Log::test() << "Testing localization: " << *local_ << std::endl;
  }

  ~LocalizationFixture<MODEL>() {}

  std::unique_ptr<const Geometry_>       resol_;
  std::unique_ptr<const oops::Variables> ctlvars_;
  std::unique_ptr<const util::DateTime>  time_;
  std::unique_ptr<Localization_>         local_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testLocalizationRandomize() {
  typedef LocalizationFixture<MODEL> Test_;
  typedef oops::Increment<MODEL>     Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());

  Test_::localization().randomize(dx);
  EXPECT(dx.norm() > 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLocalizationZero() {
  typedef LocalizationFixture<MODEL> Test_;
  typedef oops::Increment<MODEL>     Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());

  EXPECT(dx.norm() == 0.0);
  Test_::localization().multiply(dx);
  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLocalizationMultiply() {
  typedef LocalizationFixture<MODEL> Test_;
  typedef oops::Increment<MODEL>     Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx.random();

  EXPECT(dx.norm() > 0.0);
  Test_::localization().multiply(dx);
  EXPECT(dx.norm() > 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> class Localization : public oops::Test {
 public:
  Localization() {}
  virtual ~Localization() {}
 private:
  std::string testid() const override {return "test::Localization<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Localization/testLocalizationRandomize")
      { testLocalizationRandomize<MODEL>(); });
    ts.emplace_back(CASE("interface/Localization/testLocalizationZero")
      { testLocalizationZero<MODEL>(); });
    ts.emplace_back(CASE("interface/Localization/testLocalizationMultiply")
      { testLocalizationMultiply<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LOCALIZATION_H_
