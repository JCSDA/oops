/*
 * (C) Copyright 2017-2018 UCAR
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

#include "eckit/config/LocalConfiguration.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<ObsOperator_> hop(new ObsOperator_(Test_::obspace()[jj], conf[jj]));
    BOOST_CHECK(hop.get());

    hop.reset();
    BOOST_CHECK(!hop.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testSimulateObs() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    ObsOperator_ hop(Test_::obspace()[jj], conf[jj]);

    // read geovals from the file
    eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    const GeoVaLs_ gval(gconf, hop.variables());

    // initialize bias correction
    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);

    // create obsvector to hold H(x)
    ObsVector_ ovec(Test_::obspace()[jj], hop.observed());

    // call H(x), save result in the output file as @hofx
    hop.simulateObs(gval, ovec, ybias);
    ovec.save("hofx");

    const double tol = conf[jj].getDouble("tolerance");
    if (conf[jj].has("vecequiv")) {
      // if reference h(x) is saved in file as a vector, read from file
      // and compare the norm of difference to zero
      ObsVector_ ovec_ref(ovec, false);
      ovec_ref.read(conf[jj].getString("vecequiv"));
      ovec_ref -= ovec;
      const double zz = ovec_ref.rms();
      oops::Log::info() << "Vector difference between reference and computed: " <<
                           ovec_ref;
      BOOST_CHECK_SMALL(zz, tol);
    } else {
      // else compare h(x) norm to the norm from the config
      const double zz = ovec.rms();
      const double xx = conf[jj].getDouble("rmsequiv");
      BOOST_CHECK_CLOSE(xx, zz, tol);
    }
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
    ts->add(BOOST_TEST_CASE(&testSimulateObs<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSOPERATOR_H_
