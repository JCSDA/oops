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

#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/ObsAuxCovariance.h"
#include "test/TestEnvironment.h"
#include "eckit/config/Configuration.h"

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
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "ObsBiasCovariance"));
  }

  ~ObsAuxCovarianceFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration>  conf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsAuxCovarianceFixture<MODEL>   Test_;
  typedef oops::ObsAuxCovariance<MODEL>    Covariance_;

  boost::scoped_ptr<Covariance_> bias(new Covariance_(Test_::config()));
  BOOST_CHECK(bias.get());

  bias.reset();
  BOOST_CHECK(!bias.get());
}

// -----------------------------------------------------------------------------
//  void linearize(const ObsAuxControl_ &, const Geometry_ &);
//  void multiply(const ObsAuxIncrement_ &, ObsAuxIncrement_ &) const;
//  void invMult(const ObsAuxIncrement_ &, ObsAuxIncrement_ &) const;
//  void randomize(ObsAuxIncrement_ &) const;
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

template <typename MODEL> class ObsAuxCovariance : public oops::Test {
 public:
  ObsAuxCovariance() {}
  virtual ~ObsAuxCovariance() {}
 private:
  std::string testid() const {return "test::ObsAuxCovariance<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObsAuxCovariance");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXCOVARIANCE_H_
