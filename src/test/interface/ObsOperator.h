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

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

using eckit::types::is_approximately_equal;

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
    EXPECT(hop.get());

    hop.reset();
    EXPECT(!hop.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testSimulateObs() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::Locations<MODEL>         Locations_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsOperator_ hop(Test_::obspace()[jj], conf[jj]);
    eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    Locations_ locs(hop.locations(Test_::tbgn(), Test_::tend()));
    const GeoVaLs_ gval(gconf, hop.variables());

    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);

    ObsVector_ ovec(Test_::obspace()[jj], hop.observed());

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
      EXPECT(zz < tol);
    } else {
      // else compare h(x) norm to the norm from the config
      const double zz = ovec.rms();
      const double xx = conf[jj].getDouble("rmsequiv");
      EXPECT(is_approximately_equal(xx, zz, tol));
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsOperator : public oops::Test {
 public:
  ObsOperator() {}
  virtual ~ObsOperator() {}
 private:
  std::string testid() const {return "test::ObsOperator<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsOperator/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsOperator/testSimulateObs")
      { testSimulateObs<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSOPERATOR_H_
