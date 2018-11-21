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

#include <boost/scoped_ptr.hpp>

#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/ObsErrorCovariance.h"
#include "oops/interface/ObsOperator.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsTestsFixture<MODEL>  Test_;
  typedef oops::ObsErrorCovariance<MODEL>  Covar_;
  typedef oops::ObsOperator<MODEL> ObsOperator_;

  oops::instantiateObsErrorFactory<MODEL>();

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsOperator_ hop(Test_::obspace()[jj]);
    const eckit::LocalConfiguration rconf(conf[jj], "Covariance");
    boost::scoped_ptr<Covar_> R(new Covar_(rconf, Test_::obspace()[jj], hop.observed()));
    BOOST_CHECK(R.get());

    R.reset();
    BOOST_CHECK(!R.get());
  }
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
