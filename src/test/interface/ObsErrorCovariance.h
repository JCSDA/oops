/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSERRORCOVARIANCE_H_
#define TEST_INTERFACE_OBSERRORCOVARIANCE_H_

#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsErrorCovariance.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class ObsErrorCovarianceFixture : private boost::noncopyable {
  typedef oops::ObservationSpace<MODEL>        ObsSpace_;

 public:
  static const eckit::Configuration & config() {return *getInstance().conf_;}
  static const ObsSpace_ &  obspace()  {return *getInstance().obspace_;}

 private:
  static ObsErrorCovarianceFixture<MODEL>& getInstance() {
    static ObsErrorCovarianceFixture<MODEL> theObsErrorCovarianceFixture;
    return theObsErrorCovarianceFixture;
  }

  ObsErrorCovarianceFixture() {
    oops::instantiateObsErrorFactory<MODEL>();

    const util::DateTime tbgn(TestEnvironment::config().getString("window_begin"));
    const util::DateTime tend(TestEnvironment::config().getString("window_end"));

    std::vector<eckit::LocalConfiguration> obsConfs;
    TestEnvironment::config().get("Observations", obsConfs);
    BOOST_CHECK(obsConfs.size() > 0);
    const eckit::LocalConfiguration obsConf(obsConfs[0], "Observation");
    obspace_.reset(new ObsSpace_(obsConf, tbgn, tend));

    conf_.reset(new eckit::LocalConfiguration(obsConfs[0], "Covariance"));
  }

  ~ObsErrorCovarianceFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
  boost::scoped_ptr<const ObsSpace_> obspace_;
};
// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsErrorCovarianceFixture<MODEL> Test_;
  typedef oops::ObsErrorCovariance<MODEL>  Covar_;

  boost::scoped_ptr<Covar_> R(new Covar_(Test_::obspace(), Test_::config()));
  BOOST_CHECK(R.get());

  R.reset();
  BOOST_CHECK(!R.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ObsErrorCovariance : public oops::Test {
 public:
  ObsErrorCovariance() {}
  virtual ~ObsErrorCovariance() {}
 private:
  std::string testid() const {return "test::ObsErrorCovariance<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObsErrorCovariance");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSERRORCOVARIANCE_H_
