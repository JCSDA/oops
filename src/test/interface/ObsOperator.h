/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_OBSOPERATOR_H_
#define TEST_INTERFACE_OBSOPERATOR_H_

#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "test/TestEnvironment.h"
#include "test/interface/ObsTestsFixture.h"
#include "eckit/config/LocalConfiguration.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<ObsOperator_> hop(new ObsOperator_(Test_::obspace()[jj]));
    BOOST_CHECK(hop.get());

    hop.reset();
    BOOST_CHECK(!hop.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testEquiv() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsOperator_ hop(Test_::obspace()[jj]);

    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    const GeoVaLs_ gval(gconf);

    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);

    ObsVector_ ovec(Test_::obspace()[jj]);

    hop.obsEquiv(gval, ovec, ybias);

    const double zz = ovec.rms();
    const double xx = conf[jj].getDouble("rmsequiv");
    const double tol = conf[jj].getDouble("tolerance");
    BOOST_CHECK_CLOSE(xx, zz, tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ObsOperator : public oops::Test {
 public:
  ObsOperator() {}
  virtual ~ObsOperator() {}
 private:
  std::string testid() const {return "test::ObsOperator<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObsOperator");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testEquiv<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSOPERATOR_H_
