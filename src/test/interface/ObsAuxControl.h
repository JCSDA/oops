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

#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/ObsAuxControl.h"
#include "test/TestEnvironment.h"
#include "eckit/config/Configuration.h"

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
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "ObsBias"));
  }

  ~ObsAuxControlFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration>  conf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsAuxControlFixture<MODEL>   Test_;
  typedef oops::ObsAuxControl<MODEL>    ObsAux_;

  boost::scoped_ptr<ObsAux_> bias(new ObsAux_(Test_::config()));
  BOOST_CHECK(bias.get());

  bias.reset();
  BOOST_CHECK(!bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ObsAuxControlFixture<MODEL>   Test_;
  typedef oops::ObsAuxControl<MODEL>    ObsAux_;

  boost::scoped_ptr<ObsAux_> bias(new ObsAux_(Test_::config()));

  boost::scoped_ptr<ObsAux_> other(new ObsAux_(*bias));
  BOOST_CHECK(other.get());

  other.reset();
  BOOST_CHECK(!other.get());

  BOOST_CHECK(bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ObsAuxControl : public oops::Test {
 public:
  ObsAuxControl() {}
  virtual ~ObsAuxControl() {}
 private:
  std::string testid() const {return "test::ObsAuxControl<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObsAuxControl");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testCopyConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXCONTROL_H_
