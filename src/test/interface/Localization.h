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

#include <iostream>
#include <string>
#include <cmath>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "oops/runs/Test.h"
#include "oops/generic/instantiateLocalizationFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Localization.h"
#include "oops/interface/Variables.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"

namespace test {

// =============================================================================

template <typename MODEL> class LocalizationFixture : private boost::noncopyable {
  typedef oops::Localization<MODEL>   Localization_;
  typedef oops::Geometry<MODEL>       Geometry_;
  typedef oops::Variables<MODEL>      Variables_;

 public:
  static const Geometry_      & resol()        {return *getInstance().resol_;}
  static const Variables_     & ctlvars()      {return *getInstance().ctlvars_;}
  static const util::DateTime & time()         {return *getInstance().time_;}
  static const Localization_  & localization() {return *getInstance().local_;}

 private:
  static LocalizationFixture<MODEL>& getInstance() {
    static LocalizationFixture<MODEL> theLocalizationFixture;
    return theLocalizationFixture;
  }

  LocalizationFixture<MODEL>() {
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    const eckit::LocalConfiguration varConfig(TestEnvironment::config(), "Variables");
    ctlvars_.reset(new Variables_(varConfig));

    time_.reset(new util::DateTime(TestEnvironment::config().getString("TestDate")));

//  Setup the localization matrix
    oops::instantiateLocalizationFactory<MODEL>();
    const eckit::LocalConfiguration conf(TestEnvironment::config(), "Localization");
    local_.reset(new Localization_(*resol_, conf));
  }

  ~LocalizationFixture<MODEL>() {}

  boost::scoped_ptr<const Geometry_>      resol_;
  boost::scoped_ptr<const Variables_>     ctlvars_;
  boost::scoped_ptr<const util::DateTime> time_;
  boost::scoped_ptr<Localization_>        local_;
};

// =============================================================================

template <typename MODEL> void testLocalizationZero() {
  typedef LocalizationFixture<MODEL> Test_;
  typedef oops::Increment<MODEL>     Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());

  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
  Test_::localization().multiply(dx);
  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLocalizationMultiply() {
  typedef LocalizationFixture<MODEL> Test_;
  typedef oops::Increment<MODEL>     Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx.random();

  BOOST_CHECK(dx.norm() > 0.0);
  Test_::localization().multiply(dx);
  BOOST_CHECK(dx.norm() > 0.0);
}
// =============================================================================

template <typename MODEL> class Localization : public oops::Test {
 public:
  Localization() {}
  virtual ~Localization() {}
 private:
  std::string testid() const {return "test::Localization<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/Localization");

    ts->add(BOOST_TEST_CASE(&testLocalizationZero<MODEL>));
    ts->add(BOOST_TEST_CASE(&testLocalizationMultiply<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_LOCALIZATION_H_
